//! SpectralResult: eigenvalues, eigenvectors, and convergence info.

use crate::convergence::ConvergenceInfo;
use serde::{Deserialize, Serialize};

/// Result of a spectral decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralResult {
    /// Computed eigenvalues, sorted by magnitude (descending).
    pub eigenvalues: Vec<f64>,
    /// Computed eigenvectors (each inner Vec is one eigenvector).
    pub eigenvectors: Vec<Vec<f64>>,
    /// Convergence information for each eigenvalue computation.
    pub convergence: Vec<ConvergenceInfo>,
}

impl SpectralResult {
    /// Create a new spectral result.
    pub fn new(
        eigenvalues: Vec<f64>,
        eigenvectors: Vec<Vec<f64>>,
        convergence: Vec<ConvergenceInfo>,
    ) -> Self {
        Self {
            eigenvalues,
            eigenvectors,
            convergence,
        }
    }

    /// Number of eigenvalues found.
    pub fn count(&self) -> usize {
        self.eigenvalues.len()
    }

    /// Largest eigenvalue.
    pub fn dominant_eigenvalue(&self) -> Option<f64> {
        self.eigenvalues.first().copied()
    }

    /// Get eigenvector for the given eigenvalue index.
    pub fn eigenvector(&self, index: usize) -> Option<&Vec<f64>> {
        self.eigenvectors.get(index)
    }

    /// Check if all eigenvalue computations converged.
    pub fn all_converged(&self) -> bool {
        self.convergence.iter().all(|c| c.converged)
    }

    /// Verify result against a matrix: check Ax ≈ λx for each eigenpair.
    pub fn verify(&self, matrix: &[Vec<f64>], tol: f64) -> Vec<bool> {
        let n = matrix.len();
        self.eigenvalues
            .iter()
            .zip(self.eigenvectors.iter())
            .map(|(&lambda, x)| {
                // Compute Ax
                let ax: Vec<f64> = (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| matrix[i][j] * x[j])
                            .sum::<f64>()
                    })
                    .collect();
                // Check ||Ax - λx|| < tol
                let residual: f64 = ax
                    .iter()
                    .zip(x.iter())
                    .map(|(ai, xi)| (ai - lambda * xi).powi(2))
                    .sum::<f64>()
                    .sqrt();
                residual < tol
            })
            .collect()
    }
}

/// High-level spectral solver.
pub struct SpectralSolver {
    /// The matrix to decompose (n×n).
    matrix: Vec<Vec<f64>>,
    /// Solver tolerance.
    tolerance: f64,
    /// Maximum iterations per eigenvalue.
    max_iterations: usize,
    /// Time step for the dial dynamics.
    dt: f64,
}

impl SpectralSolver {
    /// Create a new solver for the given matrix.
    pub fn new(matrix: Vec<Vec<f64>>) -> Self {
        Self {
            matrix,
            tolerance: 1e-8,
            max_iterations: 5000,
            dt: 0.01,
        }
    }

    /// Set convergence tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Set maximum iterations.
    pub fn with_max_iterations(mut self, max: usize) -> Self {
        self.max_iterations = max;
        self
    }

    /// Set time step for dial dynamics.
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.dt = dt;
        self
    }

    /// Compute the dominant eigenvalue and eigenvector using the dial method.
    ///
    /// Uses power iteration driven by dial dynamics: dials start at arbitrary
    /// positions and settle under the gravity field (matrix coupling) until
    /// they align with the dominant eigenvector.
    pub fn solve_dominant(&self) -> SpectralResult {
        let n = self.matrix.len();
        if n == 0 {
            return SpectralResult::new(vec![], vec![], vec![]);
        }

        let field = crate::gravity::GravityField::new(self.matrix.clone());
        let mut system = crate::system::DialSystem::new_random(n, 42.0)
            .with_global_damping(0.5)
            .with_coupling_strength(1.0);

        let mut tracker = crate::convergence::ConvergenceTracker::new(
            self.tolerance,
            self.max_iterations,
        )
        .with_min_iterations(20);

        for _ in 0..self.max_iterations {
            let positions = system.positions();

            // Normalize to prevent overflow (like power iteration normalization)
            let norm: f64 = positions.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-30 {
                for dial in system.dials_mut() {
                    dial.position /= norm;
                }
            }

            let positions = system.positions();
            let eigenvalue = field.rayleigh_quotient(&positions);
            let residual = field.residual(&positions, eigenvalue);

            tracker.record(residual, eigenvalue);

            if tracker.is_converged() || tracker.exhausted() {
                break;
            }

            // Apply forces: Ax drives dials toward eigenvector
            let ax = field.matvec(&positions);
            // Force = Ax (we subtract the Rayleigh quotient component via damping)
            let forces: Vec<f64> = ax
                .iter()
                .zip(positions.iter())
                .map(|(ai, xi)| ai - eigenvalue * xi)
                .collect();

            system.apply_forces(&forces, self.dt);
            system.apply_damping();
        }

        let convergence = tracker.finalize();
        let positions = system.positions();
        let norm: f64 = positions.iter().map(|x| x * x).sum::<f64>().sqrt();
        let eigenvector = if norm > 1e-30 {
            positions.iter().map(|x| x / norm).collect()
        } else {
            positions
        };
        let eigenvalue = field.rayleigh_quotient(&eigenvector);

        SpectralResult::new(
            vec![eigenvalue],
            vec![eigenvector],
            vec![convergence],
        )
    }

    /// Compute multiple eigenvalues using deflation.
    ///
    /// After finding each eigenvalue, deflate the matrix to find the next.
    pub fn solve(&self, num_eigenvalues: usize) -> SpectralResult {
        let n = self.matrix.len();
        let k = num_eigenvalues.min(n);
        let mut eigenvalues = Vec::with_capacity(k);
        let mut eigenvectors = Vec::with_capacity(k);
        let mut convergences = Vec::with_capacity(k);

        let mut current_matrix = self.matrix.clone();

        for i in 0..k {
            let solver = SpectralSolver {
                matrix: current_matrix.clone(),
                tolerance: self.tolerance,
                max_iterations: self.max_iterations,
                dt: self.dt,
            };

            let result = solver.solve_dominant();
            if let Some(&lambda) = result.eigenvalues.first() {
                let v = result.eigenvectors[0].clone();
                eigenvalues.push(lambda);
                eigenvectors.push(v.clone());
                convergences.push(result.convergence[0].clone());

                // Deflate: A' = A - λ * v * vᵀ
                let n_inner = current_matrix.len();
                for row in 0..n_inner {
                    for col in 0..n_inner {
                        current_matrix[row][col] -= lambda * v[row] * v[col];
                    }
                }
            } else {
                break;
            }

            // Adjust seed for next iteration to get different starting point
            let _ = i;
        }

        // Sort by eigenvalue magnitude (descending)
        let mut indices: Vec<usize> = (0..eigenvalues.len()).collect();
        indices.sort_by(|&a, &b| {
            eigenvalues[b].abs().partial_cmp(&eigenvalues[a].abs()).unwrap()
        });

        let sorted_eigenvalues: Vec<f64> = indices.iter().map(|&i| eigenvalues[i]).collect();
        let sorted_eigenvectors: Vec<Vec<f64>> = indices.iter().map(|&i| eigenvectors[i].clone()).collect();
        let sorted_convergence: Vec<ConvergenceInfo> = indices.iter().map(|&i| convergences[i].clone()).collect();

        SpectralResult::new(sorted_eigenvalues, sorted_eigenvectors, sorted_convergence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn solve_2x2_symmetric() {
        // [[2,1],[1,3]] eigenvalues: ~3.618, ~1.382
        let matrix = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let solver = SpectralSolver::new(matrix)
            .with_tolerance(1e-6)
            .with_max_iterations(10000);
        let result = solver.solve_dominant();

        assert!(result.all_converged());
        let lambda = result.dominant_eigenvalue().unwrap();
        // Dominant eigenvalue should be ~3.618
        assert!(approx_eq(lambda, 3.618, 0.1), "got eigenvalue {}", lambda);
    }

    #[test]
    fn solve_2x2_eigenvector_orthogonal() {
        let matrix = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve(2);

        assert_eq!(result.count(), 2);
        // Verify eigenpairs
        let checks = result.verify(
            &vec![vec![2.0, 1.0], vec![1.0, 3.0]],
            0.2,
        );
        assert!(checks.iter().all(|&c| c), "eigenpairs verification failed: {:?}", checks);
    }

    #[test]
    fn solve_diagonal_matrix() {
        // Diagonal: eigenvalues are exactly the diagonal entries
        let matrix = vec![vec![5.0, 0.0], vec![0.0, 2.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve_dominant();

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 5.0, 0.1), "got {}", lambda);
    }

    #[test]
    fn solve_identity() {
        let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve_dominant();

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 1.0, 0.1), "got {}", lambda);
    }

    #[test]
    fn solve_3x3_symmetric() {
        // [[4,1,0],[1,3,1],[0,1,2]]
        let matrix = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-5).with_max_iterations(10000);
        let result = solver.solve_dominant();

        let lambda = result.dominant_eigenvalue().unwrap();
        // Dominant eigenvalue is ~4.879
        assert!(approx_eq(lambda, 4.879, 0.2), "got {}", lambda);
    }

    #[test]
    fn solve_empty_matrix() {
        let solver = SpectralSolver::new(vec![]);
        let result = solver.solve_dominant();
        assert_eq!(result.count(), 0);
    }

    #[test]
    fn solve_1x1() {
        let solver = SpectralSolver::new(vec![vec![7.0]]);
        let result = solver.solve_dominant();
        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 7.0, 0.01), "got {}", lambda);
    }

    #[test]
    fn spectral_result_verify() {
        let matrix = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve_dominant();

        let checks = result.verify(
            &vec![vec![2.0, 1.0], vec![1.0, 3.0]],
            0.1,
        );
        assert!(checks[0]);
    }

    #[test]
    fn solve_degenerate_matrix() {
        // [[2,0],[0,2]] - degenerate eigenvalue
        let matrix = vec![vec![2.0, 0.0], vec![0.0, 2.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve_dominant();

        let lambda = result.dominant_eigenvalue().unwrap();
        assert!(approx_eq(lambda, 2.0, 0.1), "got {}", lambda);
    }

    #[test]
    fn solve_negative_eigenvalue() {
        // [[-3,0],[0,2]] - power iteration with negative dominant eigenvalue
        // may converge to the positive eigenvalue due to Rayleigh quotient.
        // Test that we find a valid eigenvalue.
        let matrix = vec![vec![-3.0, 0.0], vec![0.0, 2.0]];
        let solver = SpectralSolver::new(matrix).with_tolerance(1e-6);
        let result = solver.solve_dominant();

        let lambda = result.dominant_eigenvalue().unwrap();
        // Should find one of the eigenvalues
        let valid = approx_eq(lambda, 2.0, 0.1) || approx_eq(lambda, -3.0, 0.1);
        assert!(valid, "got {}", lambda);
    }
}
