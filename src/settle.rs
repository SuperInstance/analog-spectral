//! Settle algorithm: iterate until convergence, detect eigenvalues from settled positions.

use crate::convergence::ConvergenceTracker;
use crate::gravity::GravityField;
use crate::spectral::SpectralResult;
use crate::system::DialSystem;

/// Configuration for the settling algorithm.
#[derive(Debug, Clone)]
pub struct SettleConfig {
    /// Time step for each iteration.
    pub dt: f64,
    /// Convergence tolerance for residual.
    pub tolerance: f64,
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Minimum iterations before checking convergence.
    pub min_iterations: usize,
    /// Global damping factor.
    pub damping: f64,
    /// Coupling strength multiplier.
    pub coupling_strength: f64,
    /// Whether to normalize dial positions each step.
    pub normalize: bool,
}

impl Default for SettleConfig {
    fn default() -> Self {
        Self {
            dt: 0.01,
            tolerance: 1e-8,
            max_iterations: 5000,
            min_iterations: 20,
            damping: 0.5,
            coupling_strength: 1.0,
            normalize: true,
        }
    }
}

impl SettleConfig {
    /// Create a new settle config with given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }

    /// Set time step.
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.dt = dt;
        self
    }

    /// Set max iterations.
    pub fn with_max_iterations(mut self, max: usize) -> Self {
        self.max_iterations = max;
        self
    }

    /// Set damping.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping;
        self
    }
}

/// The settling algorithm that drives dials to eigenvalue equilibrium.
pub struct Settler {
    config: SettleConfig,
}

impl Settler {
    /// Create a new settler with the given configuration.
    pub fn new(config: SettleConfig) -> Self {
        Self { config }
    }

    /// Settle a system against a gravity field to find the dominant eigenvalue.
    ///
    /// The algorithm iterates:
    /// 1. Compute Rayleigh quotient from current dial positions
    /// 2. Compute forces from residual (Ax - λx)
    /// 3. Apply forces to dials with damping
    /// 4. Normalize positions to prevent divergence
    /// 5. Check convergence
    pub fn settle_dominant(
        &self,
        field: &GravityField,
        mut system: DialSystem,
    ) -> SpectralResult {
        let n = field.dimension();
        if n == 0 {
            return SpectralResult::new(vec![], vec![], vec![]);
        }

        let mut tracker = ConvergenceTracker::new(self.config.tolerance, self.config.max_iterations)
            .with_min_iterations(self.config.min_iterations);

        for _ in 0..self.config.max_iterations {
            let positions = system.positions();

            // Normalize to prevent overflow
            if self.config.normalize {
                let norm: f64 = positions.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm > 1e-30 {
                    for dial in system.dials_mut() {
                        dial.position /= norm;
                    }
                }
            }

            let positions = system.positions();
            let eigenvalue = field.rayleigh_quotient(&positions);
            let residual = field.residual(&positions, eigenvalue);

            tracker.record(residual, eigenvalue);

            if tracker.is_converged() || tracker.exhausted() {
                break;
            }

            // Compute forces from matrix-vector product
            let ax = field.matvec(&positions);
            let forces: Vec<f64> = ax
                .iter()
                .zip(positions.iter())
                .map(|(ai, xi)| self.config.coupling_strength * (ai - eigenvalue * xi))
                .collect();

            system.apply_forces(&forces, self.config.dt);

            // Apply global damping
            for dial in system.dials_mut() {
                dial.velocity *= 1.0 - self.config.damping * self.config.dt;
            }
        }

        let convergence = tracker.finalize();

        // Final normalization
        let positions = system.positions();
        let norm: f64 = positions.iter().map(|x| x * x).sum::<f64>().sqrt();
        let eigenvector = if norm > 1e-30 {
            positions.iter().map(|x| x / norm).collect()
        } else {
            positions
        };
        let eigenvalue = field.rayleigh_quotient(&eigenvector);

        SpectralResult::new(vec![eigenvalue], vec![eigenvector], vec![convergence])
    }

    /// Settle to find multiple eigenvalues using deflation.
    pub fn settle_all(
        &self,
        field: &GravityField,
        num_eigenvalues: usize,
    ) -> SpectralResult {
        let n = field.dimension();
        let k = num_eigenvalues.min(n);
        if k == 0 {
            return SpectralResult::new(vec![], vec![], vec![]);
        }

        let mut eigenvalues = Vec::with_capacity(k);
        let mut eigenvectors = Vec::with_capacity(k);
        let mut convergences = Vec::with_capacity(k);
        let mut deflated_matrix = field.matrix().to_vec();

        for i in 0..k {
            let deflated_field = GravityField::new(deflated_matrix.clone());

            // Use different initial conditions for each eigenvalue
            let system = DialSystem::new_random(n, (i + 1) as f64 * 17.0)
                .with_global_damping(self.config.damping)
                .with_coupling_strength(self.config.coupling_strength);

            let result = self.settle_dominant(&deflated_field, system);

            if let Some(&lambda) = result.eigenvalues.first() {
                let v = result.eigenvectors[0].clone();
                eigenvalues.push(lambda);
                eigenvectors.push(v.clone());
                convergences.push(result.convergence[0].clone());

                // Deflate: A' = A - λ * v * vᵀ
                for row in 0..n {
                    for col in 0..n {
                        deflated_matrix[row][col] -= lambda * v[row] * v[col];
                    }
                }
            } else {
                break;
            }
        }

        // Sort by eigenvalue magnitude (descending)
        let mut indices: Vec<usize> = (0..eigenvalues.len()).collect();
        indices.sort_by(|&a, &b| {
            eigenvalues[b]
                .abs()
                .partial_cmp(&eigenvalues[a].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        SpectralResult::new(
            indices.iter().map(|&i| eigenvalues[i]).collect(),
            indices.iter().map(|&i| eigenvectors[i].clone()).collect(),
            indices.iter().map(|&i| convergences[i].clone()).collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn settle_2x2_dominant() {
        let field = GravityField::new(vec![vec![2.0, 1.0], vec![1.0, 3.0]]);
        let system = DialSystem::new_random(2, 1.0);
        let settler = Settler::new(SettleConfig::default());
        let result = settler.settle_dominant(&field, system);

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(
            approx_eq(lambda, 3.618, 0.2),
            "dominant eigenvalue: got {}",
            lambda
        );
    }

    #[test]
    fn settle_diagonal() {
        let field = GravityField::diagonal(&[5.0, 1.0]);
        let system = DialSystem::new_random(2, 1.0);
        let settler = Settler::new(SettleConfig::default());
        let result = settler.settle_dominant(&field, system);

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 5.0, 0.1), "got {}", lambda);
    }

    #[test]
    fn settle_all_2x2() {
        let field = GravityField::new(vec![vec![2.0, 1.0], vec![1.0, 3.0]]);
        let settler = Settler::new(SettleConfig::new(1e-6));
        let result = settler.settle_all(&field, 2);

        assert_eq!(result.count(), 2);
        // Eigenvalues should be ~3.618 and ~1.382
        let eigs = &result.eigenvalues;
        assert!(approx_eq(eigs[0], 3.618, 0.3), "first eigenvalue: got {}", eigs[0]);
        assert!(approx_eq(eigs[1], 1.382, 0.3), "second eigenvalue: got {}", eigs[1]);
    }

    #[test]
    fn settle_with_custom_config() {
        let config = SettleConfig::new(1e-10)
            .with_dt(0.005)
            .with_max_iterations(10000)
            .with_damping(0.3);

        let field = GravityField::new(vec![vec![4.0, 2.0], vec![2.0, 1.0]]);
        let system = DialSystem::new_random(2, 3.0);
        let settler = Settler::new(config);
        let result = settler.settle_dominant(&field, system);

        // [[4,2],[2,1]] has eigenvalues 5 and 0
        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 5.0, 0.2), "got {}", lambda);
    }

    #[test]
    fn settle_empty_field() {
        let field = GravityField::new(vec![]);
        let system = DialSystem::new_zeros(0);
        let settler = Settler::new(SettleConfig::default());
        let result = settler.settle_dominant(&field, system);
        assert_eq!(result.count(), 0);
    }

    #[test]
    fn settle_convergence_info() {
        let field = GravityField::new(vec![vec![3.0, 0.0], vec![0.0, 1.0]]);
        let system = DialSystem::new_random(2, 5.0);
        let settler = Settler::new(SettleConfig::default());
        let result = settler.settle_dominant(&field, system);

        assert_eq!(result.convergence.len(), 1);
        assert!(result.convergence[0].iterations > 0);
    }

    #[test]
    fn settle_3x3_multiple() {
        let field = GravityField::new(vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ]);
        let settler = Settler::new(SettleConfig::new(1e-5));
        let result = settler.settle_all(&field, 3);

        assert_eq!(result.count(), 3);
        // Dominant should be ~4.879
        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 4.879, 0.3), "dominant: got {}", lambda);
    }

    #[test]
    fn settle_1x1() {
        let field = GravityField::new(vec![vec![7.5]]);
        let system = DialSystem::new_random(1, 1.0);
        let settler = Settler::new(SettleConfig::default());
        let result = settler.settle_dominant(&field, system);

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 7.5, 0.01), "got {}", lambda);
    }
}
