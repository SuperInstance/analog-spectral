use crate::analog_dial::AnalogDial;

/// N coupled dials settle to the eigenvectors of the coupling matrix.
/// This IS Jacobi eigenvalue computation, but phrased as physical dynamics.
pub struct DialBank {
    pub dials: Vec<AnalogDial>,
    pub coupling: Vec<Vec<f64>>,
}

impl DialBank {
    pub fn new(n: usize, coupling: Vec<Vec<f64>>) -> DialBank {
        let mut dials = Vec::with_capacity(n);
        // Initialize dials at positions derived from the coupling matrix diagonal
        for i in 0..n {
            let diag = if i < coupling.len() && i < coupling[i].len() {
                coupling[i][i]
            } else {
                1.0
            };
            // Start from a small perturbation off-diagonal to seed convergence
            let mut dial = AnalogDial::new(diag, 10.0, 0.5);
            dial.position = if i == 0 { diag + 0.1 } else { diag - 0.1 * (i as f64) };
            dial.velocity = 0.0;
            dials.push(dial);
        }
        DialBank { dials, coupling }
    }

    /// All dials step simultaneously (parallel analog computation).
    pub fn step(&mut self, dt: f64) {
        let n = self.dials.len();

        // Compute forces from coupling before updating any dial
        let mut forces: Vec<f64> = vec![0.0; n];
        for i in 0..n {
            // Self-restoring force toward own setpoint
            let displacement = self.dials[i].position - self.dials[i].setpoint;
            forces[i] += -self.dials[i].gravity * displacement;

            // Coupling forces from other dials
            for j in 0..n {
                if i != j && i < self.coupling.len() && j < self.coupling[i].len() {
                    // Off-diagonal coupling pulls dial i toward weighted average with j
                    let coupling_strength = self.coupling[i][j];
                    let diff = self.dials[j].position - self.dials[i].position;
                    forces[i] += coupling_strength * diff * 0.5;
                }
            }

            // Damping
            forces[i] += -self.dials[i].friction * self.dials[i].velocity;
        }

        // Apply forces
        for i in 0..n {
            let accel = forces[i] / self.dials[i].mass;
            self.dials[i].velocity += accel * dt;
            self.dials[i].position += self.dials[i].velocity * dt;
        }
    }

    /// Settle all dials. Returns total steps.
    pub fn settle(&mut self, dt: f64, tolerance: f64) -> usize {
        let max_steps = 10_000_000;
        for i in 0..max_steps {
            let prev: Vec<f64> = self.dials.iter().map(|d| d.position).collect();
            self.step(dt);

            let max_change = self.dials.iter().zip(prev.iter())
                .map(|(d, p)| (d.position - p).abs())
                .fold(0.0f64, f64::max);

            if max_change < tolerance {
                // Also check velocities are small
                let max_vel = self.dials.iter()
                    .map(|d| d.velocity.abs())
                    .fold(0.0f64, f64::max);
                if max_vel < tolerance {
                    return i + 1;
                }
            }
        }
        max_steps
    }

    /// Read the settled positions — these ARE the eigenvector components.
    pub fn read_positions(&self) -> Vec<f64> {
        self.dials.iter().map(|d| d.position).collect()
    }

    /// Estimate eigenvalue from settled state using Rayleigh quotient.
    pub fn eigenvalue_estimate(&self) -> f64 {
        let pos = self.read_positions();
        let n = pos.len();

        // Rayleigh quotient: (x^T A x) / (x^T x)
        let mut numerator = 0.0;
        let mut denominator = 0.0;

        for i in 0..n {
            denominator += pos[i] * pos[i];
            for j in 0..n {
                if i < self.coupling.len() && j < self.coupling[i].len() {
                    numerator += pos[i] * self.coupling[i][j] * pos[j];
                }
            }
        }

        if denominator.abs() > 1e-15 {
            numerator / denominator
        } else {
            0.0
        }
    }

    /// Verify eigenvector by computing residual ||Ax - λx|| / ||x||.
    pub fn verify_eigenvector(&self, adj: &[Vec<f64>]) -> f64 {
        let pos = self.read_positions();
        let n = pos.len();
        let lambda = self.eigenvalue_estimate();

        // Compute Ax
        let mut ax = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                if i < adj.len() && j < adj[i].len() {
                    ax[i] += adj[i][j] * pos[j];
                }
            }
        }

        // Compute ||Ax - λx|| / ||x||
        let mut residual = 0.0;
        let mut norm_x = 0.0;
        for i in 0..n {
            residual += (ax[i] - lambda * pos[i]).powi(2);
            norm_x += pos[i].powi(2);
        }

        if norm_x > 1e-15 {
            residual.sqrt() / norm_x.sqrt()
        } else {
            f64::INFINITY
        }
    }
}
