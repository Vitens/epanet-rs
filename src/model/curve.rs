//! Pump head/efficiency, valve and volume curves

use crate::constants::*;
use crate::error::InputError;
use crate::model::units::{Cfs, FlowUnits, Ft, UnitSystem};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Curve {
    pub id: Box<str>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Curve {
    // Compute the intercept and slope of the curve at a given flow rate
    pub fn coefficients(&self, q: f64) -> (f64, f64) {
        // find the index of the curve segment that contains the flow
        // clamp the index to the valid range
        let x2 = self
            .x
            .partition_point(|&x| x < q)
            .max(1)
            .min(self.x.len() - 1);
        let x1 = x2 - 1;
        let y1 = self.y[x1];
        let y2 = self.y[x2];
        let r = (y2 - y1) / (self.x[x2] - self.x[x1]);
        let h0 = y1 - r * self.x[x1];
        (h0, r)
    }
}

#[derive(Debug, Clone)]
pub struct ValveCurve {
    pub curve: Curve,
}

impl ValveCurve {
    pub fn new(
        curve: &Curve,
        flow_units: &FlowUnits,
        system: &UnitSystem,
    ) -> Result<Self, InputError> {
        // clone the curve and convert the flow and head values to the standard units (CFS and Feet)
        let flows = curve
            .x
            .iter()
            .map(|x| x / flow_units.per_cfs())
            .collect::<Vec<f64>>();
        let heads = curve
            .y
            .iter()
            .map(|y| y / system.per_feet())
            .collect::<Vec<f64>>();

        let converted_curve = Curve {
            id: curve.id.clone(),
            x: flows,
            y: heads,
        };

        // validate the curve to ensure the head is decreasing and the flow is increasing monotonically
        Ok(Self {
            curve: converted_curve,
        })
    }
    pub fn coefficients(&self, q: f64) -> (f64, f64) {
        self.curve.coefficients(q)
    }
}

#[derive(Debug, Clone)]
pub struct HeadCurve {
    pub flows: Vec<Cfs>,
    pub heads: Vec<Ft>,
    pub curve_type: HeadCurveType,
    pub statistics: HeadCurveStatistics,
}

#[derive(Debug, Clone, Deserialize, Serialize, Eq, PartialEq)]
pub enum HeadCurveType {
    SinglePoint,
    ThreePointWithShutoff,
    Custom,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct HeadCurveStatistics {
    pub h_max: Ft,      // maximum head
    pub h_shutoff: Ft,  // shutoff head
    pub q_max: Cfs,     // maximum flow
    pub q_initial: Cfs, // design flow (= initial flow)
    pub r: f64,         // flow coefficient
    pub n: f64,         // pump exponent
}

impl HeadCurve {
    pub fn new(
        curve: &Curve,
        flow_units: &FlowUnits,
        system: &UnitSystem,
    ) -> Result<Self, InputError> {
        // convert the flow and head values to the standard units (CFS and Feet)
        let flows = curve.x.iter().map(|x| x / flow_units.per_cfs()).collect();
        let heads = curve.y.iter().map(|y| y / system.per_feet()).collect();

        // validate the curve to ensure the head is decreasing and the flow is increasing monotonically
        if !Self::validate_curve(&flows, &heads) {
            return Err(InputError::new(
                "Invalid head curve: Head is not decreasing or flow is not increasing monotonically",
            ));
        }

        // precompute the curve statistics
        let statistics = Self::compute_curve_statistics(&flows, &heads)?;

        let curve_type = match (flows.len(), flows[0] == 0.0) {
            (1, _) => HeadCurveType::SinglePoint,
            (3, true) => HeadCurveType::ThreePointWithShutoff,
            _ => HeadCurveType::Custom,
        };

        Ok(Self {
            flows: flows,
            heads: heads,
            curve_type: curve_type,
            statistics: statistics,
        })
    }
    pub fn coefficients(&self, q: Cfs) -> (Ft, f64) {
        // find the index of the curve segment that contains the flow
        // clamp the index to the valid range
        let x2 = self
            .flows
            .partition_point(|&x| x < q)
            .max(1)
            .min(self.flows.len() - 1);
        let x1 = x2 - 1;
        let y1 = self.heads[x1];
        let y2 = self.heads[x2];
        let r = (y2 - y1) / (self.flows[x2] - self.flows[x1]);
        let h0 = y1 - r * self.flows[x1];
        (h0, r)
    }

    // Calculate the maximum head, minimum head, and maximum flow from the curve
    pub fn compute_curve_statistics(
        flows: &Vec<Cfs>,
        heads: &Vec<Ft>,
    ) -> Result<HeadCurveStatistics, InputError> {
        if flows.len() == 1 {
            let q = flows[0];
            let h = heads[0];

            // compute the coefficients for the head curve
            let a = h * 4.0 / 3.0; // maximum head / shutoff head
            let b = (a - h) / (q * q); // flow coefficient 

            return Ok(HeadCurveStatistics {
                h_max: a,
                h_shutoff: a,
                q_max: q * 2.0,
                q_initial: q,
                r: b,
                n: 2.0,
            });
        }
        // three point curve with shutoff head at zero
        else if flows.len() == 3 && flows[0] == 0.0 {
            // try to calculate head curve statistics from three points
            // (EPANET powercurve method)
            let h0 = heads[0]; // shutoff head
            let h1 = heads[1]; // design head
            let h2 = heads[2]; // head at maximum flow
            let q1 = flows[1]; // design flow
            let q2 = flows[2]; // max flow

            let mut valid = true;

            if h0 < TINY || h0 - h1 < TINY || h1 - h2 < TINY || q1 < TINY || q2 - q1 < TINY {
                valid = false;
            }

            let h4 = h0 - h1;
            let h5 = h0 - h2;
            let c = (h5 / h4).ln() / (q2 / q1).ln();
            valid &= c > 0.0 && c <= 20.0;
            let b = -h4 / q1.powf(c);
            valid &= b < 0.0;

            if valid {
                return Ok(HeadCurveStatistics {
                    h_max: h0,
                    h_shutoff: h0,
                    q_max: q2,
                    q_initial: q1,
                    r: -b,
                    n: c,
                });
            } else {
                return Err(InputError::new(
                    "Invalid head curve: Head curve statistics could not be calculated",
                ));
            }
        } else {
            // return head curve statistics for a custom curve
            let q_max = flows[flows.len() - 1];
            let q_initial = (flows[0] + q_max) / 2.0;
            let h_max = heads[0];

            return Ok(HeadCurveStatistics {
                h_max: h_max,
                h_shutoff: h_max,
                q_max: q_max,
                q_initial: q_initial,
                r: 0.0,
                n: 1.0,
            });
        }
    }
    /// Validate the curve to ensure the head is decreasing and the flow is increasing monotonically
    pub fn validate_curve(flows: &Vec<Cfs>, heads: &Vec<Ft>) -> bool {
        for i in 1..heads.len() {
            if flows[i] <= flows[i - 1] || heads[i] >= heads[i - 1] {
                return false;
            }
        }
        true
    }
    /// Find intercept and slope of custom pump curve segment which contains speed adjusted flow
    pub fn custom_curve_coefficients(&self, q: Cfs, speed: f64) -> (f64, f64) {
        // speed adjust the flow
        let q_adjusted = q / speed;
        // compute the coefficients of the curve segment that contains the speed adjusted flow
        let (h0, r) = self.coefficients(q_adjusted);

        let mut hgrad = -r * speed;
        let mut hloss = -h0 * speed.powi(2) + hgrad * q;

        // use linear curve if gradient is too large or too small
        if hgrad > BIG_VALUE {
            hgrad = BIG_VALUE;
            hloss = -hgrad * q;
        }
        if hgrad < RQ_TOL {
            hgrad = RQ_TOL;
            hloss = -hgrad * q;
        }

        // return the gradient and head loss
        (hgrad, hloss)
    }

    // return hydraulic gradient and head loss
    pub fn curve_coefficients(&self, q: Cfs, speed: f64) -> (f64, f64) {
        // for a custom curve, find the slope and intercept of the curve segment that contains the speed adjusted flow
        if self.curve_type == HeadCurveType::Custom {
            return self.custom_curve_coefficients(q, speed);
        } else {
            // for single point and three point with shutoff curves, use the same formula as EPANET
            // shutoff head is negative to represent head gain
            // H = a - b * Q^n
            let h0 = speed.powi(2) * -self.statistics.h_shutoff;
            let mut n = self.statistics.n;
            if (self.statistics.n - 1.0) < TINY {
                n = 1.0;
            }
            let r = self.statistics.r * speed.powf(2.0 - n);

            // curve is nonlinear
            let (mut hgrad, mut hloss) = if n != 1.0 {
                // compute curve gradient
                let hgrad = n * r * q.powf(n - 1.0);
                // ... otherwise compute head loss from pump curve
                let hloss = h0 + hgrad * q / n;
                (hgrad, hloss)
            } else {
                let hgrad = r;
                let hloss = h0 + hgrad * q;
                (hgrad, hloss)
            };

            // use linear function for very small gradient
            if hgrad < RQ_TOL {
                hgrad = RQ_TOL;
                hloss = self.statistics.h_shutoff + hgrad * q / n;
            }

            (hgrad, hloss)
        }
    }
}
