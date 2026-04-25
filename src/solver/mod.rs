//! Low-level hydraulic solver and its per-step state, results and matrix assembly.

/*!
The hydraulic solver works by iteratively solving the system of equations for nodal heads and link flows until the flow errors are within the specified tolerance.

The solver works by taking a (SolverState) as input and returning a new SolverState after convergence.
The SolverState contains the flows, heads, demands, statuses, settings and resistances for the current step. When mutations are made to the network, the solverstate needs to be updated with the new properties by synchronizing the network and the solverstate.

*/

pub mod assembly;
pub mod hydraulicsolver;
pub mod matrix;
pub mod result;
pub mod state;
