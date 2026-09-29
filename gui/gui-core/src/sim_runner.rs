//! Background hydraulic simulation execution.
//!
//! `solve_hydraulics` can take real time on large or EPS networks, so it must
//! never run on the UI thread. This module mirrors the *public* stepping API
//! that `Simulation::solve_hydraulics` itself uses internally
//! (`run_hydraulics` / `next_hydraulic_timestep`) purely so incremental
//! progress can be reported over a channel; it does not reimplement any
//! solver or parsing logic — every step still goes through
//! `epanet_rs::simulation::Simulation`.

use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread;

use epanet_rs::model::network::Network;
use epanet_rs::simulation::Simulation;
use epanet_rs::solver::result::SolverResult;

#[derive(Debug)]
pub enum SimMessage {
    Progress(f32),
    Finished(Result<SolverResult, String>),
}

/// Handle to a running (or finished) background simulation.
pub struct SimHandle {
    receiver: Receiver<SimMessage>,
}

impl SimHandle {
    /// Spawn a background thread that solves `network` and reports progress.
    /// The network is a snapshot (cloned by the caller beforehand); the
    /// background thread owns it exclusively so the UI thread's copy in
    /// `AppState` can keep being edited/rendered without contention.
    pub fn spawn(network: Network, parallel: bool) -> Self {
        let (tx, rx) = channel();
        thread::spawn(move || run(network, parallel, &tx));
        SimHandle { receiver: rx }
    }

    /// Drain every message currently available without blocking. Call this
    /// once per frame from the UI event loop and turn the results into
    /// `Action::SimulationProgress` / `Action::SimulationFinished`.
    pub fn poll(&self) -> Vec<SimMessage> {
        self.receiver.try_iter().collect()
    }
}

fn run(network: Network, parallel: bool, tx: &Sender<SimMessage>) {
    let mut simulation = Simulation::new(network);

    let can_parallelize =
        parallel && !simulation.network.has_tanks() && !simulation.network.has_pressure_controls();

    if can_parallelize {
        // The parallel solver (rayon, over report steps) has no meaningful
        // incremental progress signal; report start/end only.
        let _ = tx.send(SimMessage::Progress(0.0));
        let result = simulation.solve_hydraulics(true).map_err(|e| e.to_string());
        let _ = tx.send(SimMessage::Progress(1.0));
        let _ = tx.send(SimMessage::Finished(result));
        return;
    }

    if let Err(e) = simulation.initialize_hydraulics() {
        let _ = tx.send(SimMessage::Finished(Err(e.to_string())));
        return;
    }

    let report_timestep = simulation
        .network
        .options
        .time_options
        .report_timestep
        .max(1);
    let duration = simulation.network.options.time_options.duration;
    let report_steps = duration / report_timestep + 1;

    let mut results = SolverResult::new(
        simulation.network.links.len(),
        simulation.network.nodes.len(),
        report_steps,
    );

    loop {
        let t = match simulation.run_hydraulics() {
            Ok(t) => t,
            Err(e) => {
                let _ = tx.send(SimMessage::Finished(Err(e.to_string())));
                return;
            }
        };
        if t % report_timestep == 0 {
            results.append(simulation.state.as_ref().unwrap(), t / report_timestep);
            let progress = if duration > 0 {
                (t as f32 / duration as f32).min(1.0)
            } else {
                1.0
            };
            let _ = tx.send(SimMessage::Progress(progress));
        }
        let dt = simulation.next_hydraulic_timestep();
        if dt == 0 {
            break;
        }
    }

    results.convert_units(
        &simulation.network.options.flow_units,
        &simulation.network.options.unit_system,
    );
    let _ = tx.send(SimMessage::Progress(1.0));
    let _ = tx.send(SimMessage::Finished(Ok(results)));
}
