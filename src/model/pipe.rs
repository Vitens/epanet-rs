//! `Pipe` link: resistance, headloss and coefficients for Darcy-Weisbach, Hazen-Williams and Chezy-Manning.

use crate::constants::*;
use crate::model::link::LinkStatus;
use crate::model::link::{LinkCoefficients, LinkTrait};
use crate::model::options::{HeadlossFormula, SimulationOptions};
use crate::model::units::{Cfs, Ft, UnitConversion, UnitSystem};

use fastapprox;
use serde::{Deserialize, Serialize};

// Constants used for computing Darcy-Weisbach friction factor (src: hydcoefs.c from EPANET 2.3)
const A1: f64 = 3.141_592_653_589_793_4e3; // 1000*PI
const A2: f64 = 1.570_796_326_794_896_7e3; // 500*PI
const A8: f64 = 4.618_413_198_590_667; // 5.74*(PI/4)^.9
const A9: f64 = -8.685_889_638_065_036e-1; // -2/ln(10)
const AB: f64 = 3.288_954_763_453_990_7e-3; // 5.74/(4000^.9)
const AC: f64 = -5.142_149_657_990_939e-3; // AA*AB

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Pipe {
    pub diameter: Ft,
    pub length: Ft,
    pub roughness: f64, // roughness is either in ft (Darcy-Weisbach) or unitless (Hazen-Williams, Chezy-Manning)
    pub minor_loss: f64,
    pub check_valve: bool,
    /// Headloss formula to use for the pipe
    pub headloss_formula: HeadlossFormula,
}

const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

impl LinkTrait for Pipe {
    #[inline]
    fn coefficients(
        &self,
        q: Cfs,
        r: f64,
        _setting: f64,
        status: LinkStatus,
        _: f64,
        _: f64,
    ) -> LinkCoefficients {
        if self.check_valve && q < 0.0 {
            return LinkCoefficients::new_status(1.0 / BIG_VALUE, q, LinkStatus::TempClosed);
        }
        // for closed pipes use headloss formula hloss = BIG_VALUE * q
        if status == LinkStatus::Closed || status == LinkStatus::TempClosed {
            return LinkCoefficients::simple(1.0 / BIG_VALUE, q);
        }

        if self.headloss_formula == HeadlossFormula::DarcyWeisbach {
            let (g_inv, y) = self.dw_coefficients(q, r);
            return LinkCoefficients::simple(g_inv, y);
        }

        // take the absolute value of the flow
        let q_abs = q.abs();
        // minor loss coefficient
        let ml = self.minor_loss;
        // hydraulic exponent factor
        let n = H_EXPONENT;

        // Friction head loss gradient
        let mut hgrad = n * r * q_abs.powf(n - 1.0);
        // use linear function for very small gradient
        let mut hloss = if hgrad < RQ_TOL {
            hgrad = RQ_TOL;
            hgrad * q_abs
        } else {
            hgrad * q_abs / n
        };

        // contribution of minor losses
        if ml > 0.0 {
            hloss += ml * q_abs.powi(2);
            hgrad += 2.0 * ml * q_abs;
        }
        // adjust the headloss to the sign of the flow
        hloss *= q.signum();

        // return the coefficients
        LinkCoefficients::simple(1.0 / hgrad, hloss / hgrad)
    }

    fn resistance(&self) -> f64 {
        match self.headloss_formula {
            HeadlossFormula::HazenWilliams => {
                4.727 * self.roughness.powf(-1.852) * self.diameter.powf(-4.871) * self.length
            }
            HeadlossFormula::DarcyWeisbach => {
                // D-W f factor is included in the resistance calculation
                self.length
                    / 2.0
                    / 32.2
                    / self.diameter
                    / (PI * self.diameter.powi(2) / 4.0).powi(2)
            }
            HeadlossFormula::ChezyManning => {
                panic!("Chezy Manning headloss formula not yet implemented");
            }
        }
    }

    fn update_status(
        &self,
        _: f64,
        status: LinkStatus,
        _: f64,
        _: f64,
        _: f64,
    ) -> Option<LinkStatus> {
        if status == LinkStatus::TempClosed {
            return Some(LinkStatus::Open); // reopen the pipe if it was temporarily closed
        }
        None
    }

    fn initial_flow(&self) -> f64 {
        // return the initial flow of the pipe corresponding to 1 ft/s velocity
        0.25 * std::f64::consts::PI * self.diameter.powi(2)
    }
}

impl Pipe {
    /// Calculate the coefficients for the Darcy Weisbach headloss formula
    fn dw_coefficients(&self, q: f64, r: f64) -> (f64, f64) {
        let q_abs = q.abs();
        let ml = self.minor_loss;
        let e = (self.roughness / 1000.0) / self.diameter; // relative roughness (use mf to ft)
        let s = VISCOSITY * self.diameter; // kinematic viscosity * diameter
        // Laminar flow (Re <= 2000)
        // use Hagen-Poiseuille formula
        if q_abs <= A2 * s {
            let r = 16.0 * PI * s * r;
            let hloss = q * (r + ml * q_abs);
            let hgrad = r + 2.0 * ml * q_abs;

            (1.0 / hgrad, hloss / hgrad)
        } else {
            // Turbulent flow (Re > 2000)
            let (f, dfdq) = self.dw_friction_factor(q_abs, e, s);

            let r1 = f * r + ml;
            let hloss = r1 * q_abs * q;
            let hgrad = (2.0 * r1 * q_abs) + (dfdq * r * q_abs.powi(2));

            (1.0 / hgrad, hloss / hgrad)
        }
    }

    #[inline(always)]
    // Calculate the Darcy Weisbach friction factor and its derivative
    fn dw_friction_factor(&self, q: f64, e: f64, s: f64) -> (f64, f64) {
        let w = q / s;
        // Re >= 4000, use Swamee & Jain approximation

        if w >= A1 {
            // let y1 = A8 / w.powf(0.9);
            let y1 = A8 / fastapprox::fast::pow(w as f32, 0.9) as f64;
            let y2 = e / 3.7 + y1;
            // let y3 = A9 * y2.ln();
            let y3 = A9 * fastapprox::fast::ln(y2 as f32) as f64;
            let f = 1.0 / y3.powi(2);
            let dfdq = 1.8 * f * y1 * A9 / y2 / y3 / q;

            (f, dfdq)

        // Use interpolating polynomials by E. Dunlop for transition flow (2000 < Re < 4000)
        } else {
            let y2 = e / 3.7 + AB;
            // let y3 = A9 * y2.ln();
            let y3 = A9 * fastapprox::fast::ln(y2 as f32) as f64;
            let fa = 1.0 / (y3 * y3);
            let fb = (2.0 + AC / (y2 * y3)) * fa;
            let r = w / A2;
            let x1 = 7.0 * fa - fb;
            let x2 = 0.128 - 17.0 * fa + 2.5 * fb;
            let x3 = -0.128 + 13.0 * fa - (fb + fb);
            let x4 = 0.032 - 3.0 * fa + 0.5 * fb;
            let f = x1 + r * (x2 + r * (x3 + r * x4));
            let dfdq = (x2 + r * (2.0 * x3 + r * 3.0 * x4)) / s / A2;

            (f, dfdq)
        }
    }
}

impl UnitConversion for Pipe {
    fn convert_to_standard(&mut self, options: &SimulationOptions) {
        if options.unit_system == UnitSystem::US {
            self.diameter /= 12.0; // convert in to ft
        } else {
            self.diameter /= 1000.0; // convert mm to m
        }

        // convert diameter from the given unit system to feet
        self.diameter /= options.unit_system.per_feet();

        // convert length from the given unit system to feet
        self.length /= options.unit_system.per_feet();

        // if the headloss formula is Darcy Weisbach, convert roughness from the given unit system to feet
        if self.headloss_formula == HeadlossFormula::DarcyWeisbach {
            self.roughness /= options.unit_system.per_feet();
        }
    }

    fn convert_from_standard(&mut self, options: &SimulationOptions) {
        if options.unit_system == UnitSystem::US {
            self.diameter *= 12.0; // convert ft to in
        } else {
            self.diameter *= options.unit_system.per_feet() * 1000.0; // convert ft to mm
        }

        // convert length from feet to the given unit system
        self.length *= options.unit_system.per_feet();

        // if the headloss formula is Darcy Weisbach, convert roughness from feet to the given unit system
        if self.headloss_formula == HeadlossFormula::DarcyWeisbach {
            self.roughness *= options.unit_system.per_feet();
        }
    }
}
