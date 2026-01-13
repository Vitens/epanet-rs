#[derive(Debug, Clone)]
pub struct Curve {
  pub id: Box<str>,
  pub x: Vec<f64>,
  pub y: Vec<f64>,
}