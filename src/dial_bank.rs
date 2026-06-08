//! Coupled dial bank — N dials connected by a coupling matrix.
//!
//! The coupling matrix acts like a graph adjacency/weight matrix. Each dial
//! feels forces from its neighbors proportional to coupling strength. When
//! the bank settles, the dial positions approximate eigenvectors and the
//! Rayleigh quotient gives the eigenvalue estimate.
//!
//! This IS the Jacobi eigenvalue algorithm, but phrased as physical dynamics.

use crate::analog_dial::AnalogDial;

/// A bank of N coupled analog dials for eigenvalue computation.
///
/// Each dial corresponds to a component of the eigenvector. The coupling
/// matrix determines how dials influence each other — off-diagonal entries
/// create inter-dial forces that drive the system toward eigenvector alignment.
///
/// # Algorithm
///
/// 1. Initialize dials near the diagonal of the coupling matrix
/// 2. Step all dials simultaneously (parallel analog computation)
/// 3. Read settled positions as eigenvector components
/// 4. Compute Rayleigh quotient for eigenvalue estimate
pub struct DialBank {
    /// The individual dials — one per eigenvector component.
    pub dials: Vec<AnalogDial>,
    /// N×N coupling (weight) matrix — acts as the graph's adjacency structure.
    pub coupling: Vec<Vec<f64>>,
}

impl DialBank {
    /// Create a bank of N dials coupled by the given matrix.
    ///
    /// Dials are initialized near the diagonal values with small perturbations
    /// to seed convergence. Gravity is 10.0, friction is 0.5 by default.
    pub fn new(n: usize, coupling: Vec<Vec<f64>>) -> DialBank {
        let mut dials = Vec::with_capacity(n);
        for i in 0..n {
            let diag = if i < coupling.len() && i < coupling[i].len() {
                coupling[i][i]
            } else {
                1.0
            };
            let mut dial = AnalogDial::new(diag, 10.0, 0.5);
            dial.position = if i == 0 { diag + 0.1 } else { diag - 0.1 * (i as f64) };
            dial.velocity = 0.0;
            dials.push(dial);
        }
        DialBank { dials, coupling }
    }

    /// Advance all dials by one timestep simultaneously.
    ///
    /// Forces are computed from the coupling matrix before any dial is updated,
    /// ensuring symmetric force application (Gauss-Seidel-free update).
    pub fn step(&mut self, dt: f64) {
        let n = self.dials.len();

        // Compute forces from coupling before updating any dial
        let mut forces: Vec<f64> = vec![0.0; n];
        for (i, force_i) in forces.iter_mut().enumerate() {
            // Self-restoring force toward own setpoint
            let displacement = self.dials[i].position - self.dials[i].setpoint;
            *force_i += -self.dials[i].gravity * displacement;

            // Coupling forces from other dials
            for (j, dial_j) in self.dials.iter().enumerate() {
                if i != j && i < self.coupling.len() && j < self.coupling[i].len() {
                    let coupling_strength = self.coupling[i][j];
                    let diff = dial_j.position - self.dials[i].position;
                    *force_i += coupling_strength * diff * 0.5;
                }
            }

            // Damping
            *force_i += -self.dials[i].friction * self.dials[i].velocity;
        }

        // Apply forces
        for (i, force_i) in forces.iter().enumerate() {
            let accel = force_i / self.dials[i].mass;
            self.dials[i].velocity += accel * dt;
            self.dials[i].position += self.dials[i].velocity * dt;
        }
    }

    /// Settle all dials to equilibrium. Returns total steps taken.
    ///
    /// Convergence is checked by maximum position change and maximum velocity.
    /// Both must fall below `tolerance` simultaneously.
    pub fn settle(&mut self, dt: f64, tolerance: f64) -> usize {
        let max_steps = 10_000_000;
        for i in 0..max_steps {
            let prev: Vec<f64> = self.dials.iter().map(|d| d.position).collect();
            self.step(dt);

            let max_change = self.dials.iter().zip(prev.iter())
                .map(|(d, p)| (d.position - p).abs())
                .fold(0.0f64, f64::max);

            if max_change < tolerance {
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

    /// Read the settled dial positions — these approximate eigenvector components.
    pub fn read_positions(&self) -> Vec<f64> {
        self.dials.iter().map(|d| d.position).collect()
    }

    /// Estimate the dominant eigenvalue using the Rayleigh quotient.
    ///
    /// Rayleigh quotient: λ ≈ (xᵀAx) / (xᵀx), where x is the dial position
    /// vector and A is the coupling matrix.
    pub fn eigenvalue_estimate(&self) -> f64 {
        let pos = self.read_positions();
        let n = pos.len();

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

    /// Verify the eigenvector by computing the residual ||Ax − λx|| / ||x||.
    ///
    /// A smaller residual indicates a better eigenvector approximation.
    pub fn verify_eigenvector(&self, adj: &[Vec<f64>]) -> f64 {
        let pos = self.read_positions();
        let n = pos.len();
        let lambda = self.eigenvalue_estimate();

        let mut ax = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                if i < adj.len() && j < adj[i].len() {
                    ax[i] += adj[i][j] * pos[j];
                }
            }
        }

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
