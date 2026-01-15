use crate::model::link::LinkTrait;
use crate::model::options::HeadlossFormula;
use crate::model::units::{FlowUnits, UnitSystem, UnitConversion};
use crate::constants::*;


// Constants used for computing Darcy-Weisbach friction factor (src: hydcoefs.c from EPANET 2.3)
const A1 : f64 =  3.14159265358979323850e+03;   // 1000*PI
const A2 : f64 =  1.57079632679489661930e+03;   // 500*PI
// const A3 : f64 =  5.02654824574366918160e+01;   // 16*PI
// const A4 : f64 =  6.28318530717958647700e+00;   // 2*PI
const A8 : f64 =  4.61841319859066668690e+00;   // 5.74*(PI/4)^.9
const A9 : f64 = -8.68588963806503655300e-01;  // -2/ln(10)
// const AA : f64 = -1.5634601348517065795e+00;   // -2*.9*2/ln(10)
const AB : f64 =  3.28895476345399058690e-03;   // 5.74/(4000^.9)
const AC : f64 = -5.14214965799093883760e-03;  // AA*AB

#[derive(Debug, Eq, PartialEq)]
pub enum PipeStatus {
  Open,
  Closed,
  CheckValve
}
impl PipeStatus {
  pub fn from_str(status: &str) -> PipeStatus {
    match status.to_uppercase().as_str() {
      "OPEN" => PipeStatus::Open,
      "CLOSED" => PipeStatus::Closed,
      "CV" => PipeStatus::CheckValve,
      _ => panic!("Invalid pipe status")
    }
  }
}

pub struct Pipe {
  pub diameter: f64,
  pub length: f64,
  pub roughness: f64,
  pub minor_loss: f64,
  pub status: PipeStatus,
  /// Headloss formula to use for the pipe
  pub headloss_formula: HeadlossFormula
}

const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

impl LinkTrait for Pipe {
  fn coefficients(&self, q: f64, r: f64) -> (f64, f64) {

    // for closed pipes use headloss formula hloss = BIG_VALUE * q
    if self.status == PipeStatus::Closed {
      return (1.0 / BIG_VALUE, q);
    }

    if self.headloss_formula == HeadlossFormula::DarcyWeisbach {
      return self.dw_coefficients(q, r);
    }

    // take the absolute value of the flow
    let q_abs = q.abs();
    // minor loss coefficient
    let ml = self.minor_loss;
    // hydraulic exponent factor
    let n = H_EXPONENT;

    // Friction head loss gradient
    let mut hgrad = n * r * q_abs.powf(n-1.0);
    // Headloss
    let mut hloss = hgrad * q_abs / n;

    // contribution of minor losses
    if ml > 0.0 {
      hloss += ml * q_abs.powi(2);
      hgrad += 2.0 * ml * q_abs;
    }
    // adjust the headloss to the sign of the flow
    hloss *= q.signum();

    // return the coefficients
    (1.0 / hgrad, hloss/hgrad)
  }

  fn resistance(&self) -> f64 {
    match self.headloss_formula {
      HeadlossFormula::HazenWilliams => {
        4.727 * self.roughness.powf(-1.852) * self.diameter.powf(-4.87) * self.length
      }
      HeadlossFormula::DarcyWeisbach => {
        // D-W f factor is included in the resistance calculation
        self.length / 2.0 / 32.2 / self.diameter / (PI * self.diameter.powi(2) / 4.0).powi(2)
      }
      HeadlossFormula::ChezyManning => {
        panic!("Chezy Manning headloss formula not yet implemented");
      }
    }
  }
}

impl Pipe {
  /// Calculate the coefficients for the Darcy Weisbach headloss formula
  fn dw_coefficients(&self, q: f64, r: f64) -> (f64, f64) {

    let q_abs = q.abs();
    let ml = self.minor_loss;
    let e = (self.roughness / 1000.0) / self.diameter; // relative roughness (use mf to ft)
    let s = VISCOSITY * self.diameter;      // kinematic viscosity * diameter

    // Laminar flow (Re <= 2000)
    // use Hagen-Poiseuille formula
    if q_abs <= A2 * s {
      let r = 16.0 * PI * s * r;
      let hloss = q * (r + ml * q_abs);
      let hgrad = r + 2.0 * ml * q_abs;

      (1.0/hgrad, hloss/hgrad)

    } else {
      // Turbulent flow (Re > 2000)
      let (f, dfdq) = self.dw_friction_factor(q_abs, e, s);

      let r1 = f * r + ml;
      let hloss = r1 * q_abs * q;
      let hgrad = (2.0 * r1 * q_abs) + (dfdq * r * q_abs.powi(2));

      // println!("hloss: {}, hgrad: {}", hloss, hgrad);

      (1.0/hgrad, hloss/hgrad)

    }

  }

  #[inline(always)]
  // Calculate the Darcy Weisbach friction factor and its derivative
  fn dw_friction_factor(&self, q: f64, e: f64, s: f64) -> (f64, f64) {

    let w = q / s;
    // Re >= 4000, use Swamee & Jain approximation

    if w >= A1 {
      let y1 = A8 / w.powf(0.9);
      let y2 = e / 3.7 + y1;
      let y3 = A9 * y2.ln();
      let f = 1.0 / y3.powi(2);
      let dfdq = 1.8 * f * y1 * A9 / y2 / y3 / q;

      (f, dfdq)
    
    // Use interpolating polynomials by E. Dunlop for transition flow (2000 < Re < 4000)
    } else {
      let y2 = e / 3.7 + AB;
      let y3 = A9 * y2.ln();
      let fa = 1.0 / (y3*y3);
      let fb = (2.0 + AC / (y2*y3)) * fa;
      let r = w / A2;
      let x1 = 7.0 * fa - fb;
      let x2 = 0.128 - 17.0 * fa + 2.5 * fb;
      let x3 = -0.128 + 13.0 * fa - (fb + fb);
      let x4 = 0.032 - 3.0 * fa + 0.5 *fb;
      let f = x1 + r * (x2 + r * (x3 + r * x4));
      let dfdq = (x2 + r * (2.0 * x3 + r * 3.0 * x4)) / s / A2;

      (f, dfdq)
    }

  }
}

impl UnitConversion for Pipe {
  fn convert_units(&mut self, _flow: &FlowUnits, system: &UnitSystem, reverse: bool) {
    if system == &UnitSystem::SI {

      if reverse {
        self.diameter = self.diameter * MperFT * 1e3; // convert in to ft to mm
        self.length = self.length * MperFT; // convert ft to m
        self.roughness = self.roughness * MperFT; // convert mmft to ft
      }
      else {
        self.diameter = self.diameter / 1e3 / MperFT; // convert mm to in
        self.length = self.length / MperFT;
        self.roughness = self.roughness / MperFT; // convert mm to mmft
      }
    }
    else {
      self.diameter = self.diameter / 12.0; // convert in to ft
    }
  }
}