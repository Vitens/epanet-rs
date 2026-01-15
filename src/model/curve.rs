#[derive(Debug, Clone)]
pub struct Curve {
  pub id: Box<str>,
  pub x: Vec<f64>,
  pub y: Vec<f64>,
}


#[derive(Debug, Clone)]
pub struct HeadCurveStatistics {
  pub h_max: f64,           // maximum head
  pub h_shutoff: f64,       // shutoff head
  pub q_max: f64,           // maximum flow
  pub q_initial: f64,       // design flow (= initial flow)
  pub r: f64,               // flow coefficient
  pub n: f64,               // pump exponent
}

impl Curve {
  // Calculate the maximum head, minimum head, and maximum flow from the curve
  pub fn head_curve_statistics(&self) -> HeadCurveStatistics {

    if self.x.len() == 1 {

      let q = self.x[0];
      let h = self.y[0];

      // compute the coefficients for the head curve
      let a = h * 4.0 / 3.0; // maximum head / shutoff head
      let b = (a-h)/(q*q);  // flow coefficient 

      return HeadCurveStatistics {
        h_max: a,
        h_shutoff: a,
        q_max: q * 2.0,
        q_initial: q,
        r: b,
        n: 2.0,
      }
    }

    else {
      panic!("Only single point curves are supported for now");
    }
  }
}
