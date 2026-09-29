//! EPANET network editor — eframe/egui desktop binary.
//!
//! Usage: `epanet-gui [path/to/network.inp]`

mod app;
mod canvas;
mod command_palette;
mod panels;
mod shortcuts;

use epanet_rs::model::network::Network;
use epanet_rs::model::options::HeadlossFormula;
use epanet_rs::model::units::FlowUnits;
use gui_core::AppState;

fn main() -> eframe::Result<()> {
    let network = std::env::args()
        .nth(1)
        .and_then(|path| match Network::from_file(&path) {
            Ok(network) => Some(network),
            Err(e) => {
                eprintln!("Failed to load '{path}': {e}");
                None
            }
        })
        .unwrap_or_else(|| Network::new(FlowUnits::LPS, HeadlossFormula::DarcyWeisbach));

    let mut state = AppState::new(network);
    if let Some(bounds) = state.content_bounds() {
        state.camera.fit_to_bounds(bounds);
    }

    let native_options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_title("EPANET Network Editor"),
        ..Default::default()
    };

    eframe::run_native(
        "EPANET Network Editor",
        native_options,
        Box::new(|_cc| Ok(Box::new(app::EditorApp::new(state)))),
    )
}
