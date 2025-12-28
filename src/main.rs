use std::time::Instant;

struct Pipes {
  flow: Vec<f64>,
  diameter: Vec<f64>,
  length: Vec<f64>,
  roughness: Vec<f64>,
  headloss: Vec<f64>,
  gradient: Vec<f64>
}

impl Pipes {

  fn fast_pow_09(re: f64) -> f64 {
    let b = 0.9;
    let mut u = re.to_bits();
    // Magic constant for f64: (1 - b) * 1023 in bits
    u = (b * (u as f64 - 4606045658400000000.0) + 4606045658400000000.0) as u64;
    f64::from_bits(u)
  }

  fn headloss(&mut self) {
    for i in 0..self.flow.len() {
      let flow = self.flow[i].abs();
      let diameter = self.diameter[i];
      let length = self.length[i];
      let roughness = self.roughness[i];

      const A1: f64 = 0.7853981634; // pi/4
      const A2: f64 = -0.8685889638065037; // -2/ln(10)

      let diameter_m = diameter/1000.0;
      let area = A1 * diameter_m.powi(2);
      let velocity = flow / area / 3600.0; // velocity in m/s
      let re = velocity * diameter_m / 1.0e-6; // Reynolds number
      let roughness_m = roughness/1000.0; // roughness in m
      let epsilon_d = roughness_m / diameter_m;

      const RHO: f64 = 1000.0; // density of water in kg/m^3
      let y1 = 5.74 / Self::fast_pow_09(re);
      let y2 = epsilon_d / 3.7 + y1;
      let y3 = A2 * y2.ln(); // todo: improve this
      let f = 1.0 / (y3 * y3);
      let df_dq = 1.8 * f * y1 * A2 / y2 / y3 / flow;
      let headloss = f * length * velocity * velocity * RHO / (2.0 * diameter_m) / 9810.0;
      self.headloss[i] = headloss;
      self.gradient[i] = (2.0 * headloss / flow) + (headloss / f * df_dq);
    }
  }

  }

fn main() {

  let iterations = 10_000_000;
  println!("Creating {} pipes", iterations);
  let start = Instant::now();
  let mut pipes = Pipes {
    flow: vec![0.0; iterations],
    diameter: vec![300.0; iterations],
    length: vec![100.0; iterations],
    roughness: vec![0.25; iterations],
    headloss: vec![0.0; iterations],
    gradient: vec![0.0; iterations]
  };
  for i in 0..iterations {
    pipes.flow[i] = i as f64 + 1.0;
  }
  let duration = start.elapsed();
  println!("Time taken for creation: {:?}", duration);
  let start = Instant::now();

  pipes.headloss();
  let tot = pipes.headloss.iter().sum::<f64>();
  let duration = start.elapsed();
  println!("Time taken for headloss: {:?}", duration);
  println!("Iterations per second: {}", iterations as f64 / duration.as_secs_f64());
  println!("Total headloss: {}", tot);
}
