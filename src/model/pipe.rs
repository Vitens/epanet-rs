use crate::model::link::LinkTrait;

pub enum PipeStatus {
  Open,
  Closed,
  CheckValve
}

pub struct Pipe {
  pub diameter: f64,
  pub length: f64,
  pub roughness: f64,
  pub minor_loss: f64,
  pub status: PipeStatus,
}

const H_EXPONENT: f64 = 1.852; // Hazen-Williams exponent

impl LinkTrait for Pipe {
  fn coefficients(&self, q: f64, resistance: f64) -> (f64, f64) {
    let q_abs = q.abs().max(1e-8);
    let r = resistance;
    let m = self.minor_loss;
    let n = H_EXPONENT;

    let q_pow = q_abs.powf(n - 1.0);

    let g = n * r * q_pow + 2.0 * m * q_abs;
    let g_inv = 1.0 / g;

    let y = (r * q_abs * q_pow + m * q_abs.powi(2)) * q.signum();

    (g_inv, y)
  }
  fn resistance(&self) -> f64 {
    0.0
  }
}