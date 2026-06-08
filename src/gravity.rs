//! Gravity field: coupling forces between dials.
//!
//! The gravity field implements the matrix-vector product that drives dial
//! dynamics. Equilibrium positions correspond to eigenvectors.

use crate::system::DialSystem;

/// Represents the coupling forces in the analog spectral system.
///
/// The "gravity" field encodes the matrix A as coupling forces: when dials
/// are at positions forming a vector x, the force on each dial is A*x - x,
/// which drives the system toward eigenvector equilibria.
#[derive(Debug, Clone)]
pub struct GravityField {
    /// The matrix whose eigenvalues we seek (row-major, n×n).
    matrix: Vec<Vec<f64>>,
    /// Matrix dimension.
    n: usize,
}

impl GravityField {
    /// Create a gravity field from an n×n matrix.
    pub fn new(matrix: Vec<Vec<f64>>) -> Self {
        let n = matrix.len();
        Self { matrix, n }
    }

    /// Create a gravity field from a flat row-major slice.
    pub fn from_flat(data: &[f64], n: usize) -> Self {
        let matrix = data.chunks(n).map(|row| row.to_vec()).collect();
        Self { matrix, n }
    }

    /// Matrix dimension.
    pub fn dimension(&self) -> usize {
        self.n
    }

    /// Compute matrix-vector product y = A*x.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.n, "Vector dimension mismatch");
        let mut y = vec![0.0; self.n];
        for (i, yi) in y.iter_mut().enumerate().take(self.n) {
            for (j, xj) in x.iter().enumerate().take(self.n) {
                *yi += self.matrix[i][j] * xj;
            }
        }
        y
    }

    /// Compute the Rayleigh quotient λ = (xᵀAx) / (xᵀx).
    ///
    /// This gives the eigenvalue estimate for a given dial configuration.
    pub fn rayleigh_quotient(&self, x: &[f64]) -> f64 {
        let ax = self.matvec(x);
        let xtax: f64 = x.iter().zip(ax.iter()).map(|(xi, ai)| xi * ai).sum();
        let xtx: f64 = x.iter().map(|xi| xi * xi).sum();
        if xtx < 1e-30 {
            0.0
        } else {
            xtax / xtx
        }
    }

    /// Compute the residual ||Ax - λx|| for a given eigenvalue estimate λ.
    pub fn residual(&self, x: &[f64], lambda: f64) -> f64 {
        let ax = self.matvec(x);
        let res: f64 = ax
            .iter()
            .zip(x.iter())
            .map(|(ai, xi)| (ai - lambda * xi).powi(2))
            .sum();
        res.sqrt()
    }

    /// Compute forces on each dial from the gravity field.
    ///
    /// Force on dial i: F_i = (Ax)_i - λ * x_i
    /// This drives the system toward eigenvectors where F → 0.
    pub fn compute_forces(&self, system: &DialSystem, lambda: f64) -> Vec<f64> {
        let positions: Vec<f64> = system.dials().iter().map(|d| d.position).collect();
        let ax = self.matvec(&positions);
        ax.iter()
            .zip(positions.iter())
            .map(|(ai, xi)| ai - lambda * xi)
            .collect()
    }

    /// Get a reference to the underlying matrix.
    pub fn matrix(&self) -> &[Vec<f64>] {
        &self.matrix
    }

    /// Compute the trace of the matrix.
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.matrix[i][i]).sum()
    }

    /// Compute the Frobenius norm of the matrix.
    pub fn frobenius_norm(&self) -> f64 {
        let sum_sq: f64 = self
            .matrix
            .iter()
            .flat_map(|row| row.iter())
            .map(|v| v * v)
            .sum();
        sum_sq.sqrt()
    }

    /// Check if the matrix is symmetric.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                if (self.matrix[i][j] - self.matrix[j][i]).abs() > tol {
                    return false;
                }
            }
        }
        true
    }

    /// Create an n×n identity matrix field.
    pub fn identity(n: usize) -> Self {
        let matrix = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();
        Self { matrix, n }
    }

    /// Create a diagonal matrix field.
    pub fn diagonal(values: &[f64]) -> Self {
        let n = values.len();
        let matrix = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { values[i] } else { 0.0 })
                    .collect()
            })
            .collect();
        Self { matrix, n }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn matvec_identity() {
        let g = GravityField::identity(3);
        let x = vec![1.0, 2.0, 3.0];
        let y = g.matvec(&x);
        assert_eq!(y, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn matvec_2x2() {
        let g = GravityField::new(vec![vec![2.0, 1.0], vec![1.0, 3.0]]);
        let y = g.matvec(&vec![1.0, 0.0]);
        assert_eq!(y, vec![2.0, 1.0]);
    }

    #[test]
    fn rayleigh_quotient_eigenvector() {
        // [[2,0],[0,5]] has eigenvector [0,1] with eigenvalue 5
        let g = GravityField::new(vec![vec![2.0, 0.0], vec![0.0, 5.0]]);
        let x = vec![0.0, 1.0];
        let rq = g.rayleigh_quotient(&x);
        assert!(approx_eq(rq, 5.0, 1e-10));
    }

    #[test]
    fn rayleigh_quotient_zero_vector() {
        let g = GravityField::identity(2);
        assert_eq!(g.rayleigh_quotient(&[0.0, 0.0]), 0.0);
    }

    #[test]
    fn residual_at_eigenvector() {
        // [[2,0],[0,5]] with eigenvector [0,1] and eigenvalue 5
        let g = GravityField::new(vec![vec![2.0, 0.0], vec![0.0, 5.0]]);
        let x = vec![0.0, 1.0];
        let res = g.residual(&x, 5.0);
        assert!(res < 1e-10);
    }

    #[test]
    fn trace_computation() {
        let g = GravityField::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert!(approx_eq(g.trace(), 5.0, 1e-10));
    }

    #[test]
    fn frobenius_norm() {
        let g = GravityField::new(vec![vec![3.0, 4.0]]);
        // sqrt(9+16) = 5
        assert!(approx_eq(g.frobenius_norm(), 5.0, 1e-10));
    }

    #[test]
    fn is_symmetric() {
        let g = GravityField::new(vec![vec![1.0, 2.0], vec![2.0, 3.0]]);
        assert!(g.is_symmetric(1e-10));
    }

    #[test]
    fn not_symmetric() {
        let g = GravityField::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert!(!g.is_symmetric(1e-10));
    }

    #[test]
    fn from_flat() {
        let g = GravityField::from_flat(&[1.0, 0.0, 0.0, 1.0], 2);
        let y = g.matvec(&[3.0, 4.0]);
        assert_eq!(y, vec![3.0, 4.0]);
    }

    #[test]
    fn diagonal_matrix() {
        let g = GravityField::diagonal(&[2.0, 5.0]);
        let y = g.matvec(&[1.0, 1.0]);
        assert_eq!(y, vec![2.0, 5.0]);
    }
}
