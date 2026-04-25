//! `Junction` node: a demand point with optional emitter and pressure-dependent demand coefficients.

use crate::constants::{BIG_VALUE, MperFT, PSIperFT, RQ_TOL, SMALL_VALUE};
use crate::model::options::SimulationOptions;
use crate::model::units::{Cfs, UnitConversion, UnitSystem};

use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Junction {
    pub basedemand: Cfs,
    pub emitter_coefficient: f64,
    pub pattern: Option<Box<str>>,
    #[serde(skip)]
    pub pattern_index: Option<usize>,
}

impl Junction {
    /// Compute the coefficients for the demand function
    /// demand: current demand of junction
    /// demand_full: full demand
    /// dp: pressure range for demand function
    /// n: 1 / pressure exponent
    pub fn demand_coefficients(
        &self,
        demand: Cfs,
        demand_full: Cfs,
        dp: f64,
        n: f64,
    ) -> (f64, f64) {
        // calculate factor
        let factor = demand / demand_full;

        // use a large value for the gradient and headloss to prevent negative demand
        if factor <= 0.0 {
            return (1.0 / BIG_VALUE, BIG_VALUE * demand);
        } else if factor < 1.0 {
            let hgrad = n * dp * factor.powf(n - 1.0) / demand_full;
            // use a linear function for very small gradient
            if hgrad < RQ_TOL {
                return (1.0 / RQ_TOL, demand * RQ_TOL);
            }

            let y = hgrad * demand / n;
            return (1.0 / hgrad, y);
        } else {
            // use a large value for the gradient and headloss to prevent demand above full value
            return (1.0 / BIG_VALUE, dp + BIG_VALUE * (demand - demand_full));
        }
    }
    /// Compute the coefficients for the emitter function
    /// q: flow through the emitter
    /// q_exp: emitter exponent
    pub fn emitter_coefficients(&self, q: Cfs, q_exp: f64) -> (f64, f64) {
        // get the emitter coefficient
        let ke = self.emitter_coefficient.max(SMALL_VALUE);

        // compute the gradient of headloss through emitter
        let hgrad = q_exp * ke * q.abs().powf(q_exp - 1.0);

        if hgrad < RQ_TOL {
            let hgrad = RQ_TOL;
            return (1.0 / hgrad, q * hgrad);
        }

        let y = hgrad * q / q_exp;

        return (1.0 / hgrad, y);
    }
}

impl UnitConversion for Junction {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        self.basedemand = self.basedemand / options.flow_units.per_cfs();
        // convert the emitter coefficient to the standard unit system
        if self.emitter_coefficient > 0.0 {
            let pressure_factor = if options.unit_system == UnitSystem::US {
                PSIperFT
            } else {
                MperFT
            };
            let conversion_factor =
                options.flow_units.per_cfs().powf(options.emitter_exponent) / pressure_factor;

            self.emitter_coefficient =
                conversion_factor / self.emitter_coefficient.powf(options.emitter_exponent);
        }
    }
    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        self.basedemand = self.basedemand * options.flow_units.per_cfs();
        // convert the emitter coefficient from the standard unit system
        if self.emitter_coefficient > 0.0 {
            let pressure_factor = if options.unit_system == UnitSystem::US {
                PSIperFT
            } else {
                MperFT
            };
            let conversion_factor =
                options.flow_units.per_cfs().powf(options.emitter_exponent) / pressure_factor;

            self.emitter_coefficient =
                (conversion_factor / self.emitter_coefficient).powf(1.0 / options.emitter_exponent);
        }
    }
}
