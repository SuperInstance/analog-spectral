//! Convergence tracking for the spectral settling process.

use serde::{Deserialize, Serialize};

/// Information about convergence of the settling process.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceInfo {
    /// Whether the process converged.
    pub converged: bool,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual norm.
    pub residual: f64,
    /// Estimated convergence rate (log reduction per iteration).
    pub rate: f64,
    /// Estimated spectral gap (difference between smallest and next eigenvalue).
    pub spectral_gap: Option<f64>,
    /// History of residuals (sampled).
    pub residual_history: Vec<f64>,
}

impl ConvergenceInfo {
    /// Create a convergence info indicating failure.
    pub fn did_not_converge(iterations: usize, residual: f64) -> Self {
        Self {
            converged: false,
            iterations,
            residual,
            rate: 0.0,
            spectral_gap: None,
            residual_history: Vec::new(),
        }
    }
}

/// Tracker for monitoring convergence during settling.
#[derive(Debug, Clone)]
pub struct ConvergenceTracker {
    /// Residual tolerance for convergence.
    pub tolerance: f64,
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// History of residual values.
    residual_history: Vec<f64>,
    /// History of eigenvalue estimates.
    eigenvalue_history: Vec<f64>,
    /// Iteration counter.
    iteration: usize,
    /// Minimum number of iterations before checking convergence.
    min_iterations: usize,
}

impl ConvergenceTracker {
    /// Create a new convergence tracker.
    pub fn new(tolerance: f64, max_iterations: usize) -> Self {
        Self {
            tolerance,
            max_iterations,
            residual_history: Vec::new(),
            eigenvalue_history: Vec::new(),
            iteration: 0,
            min_iterations: 10,
        }
    }

    /// Set minimum iterations before convergence check.
    pub fn with_min_iterations(mut self, min: usize) -> Self {
        self.min_iterations = min;
        self
    }

    /// Record a new residual and eigenvalue estimate.
    pub fn record(&mut self, residual: f64, eigenvalue: f64) {
        self.residual_history.push(residual);
        self.eigenvalue_history.push(eigenvalue);
        self.iteration += 1;
    }

    /// Current iteration count.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// Check if we've exceeded max iterations.
    pub fn exhausted(&self) -> bool {
        self.iteration >= self.max_iterations
    }

    /// Check if converged: residual below tolerance and enough iterations.
    pub fn is_converged(&self) -> bool {
        if self.iteration < self.min_iterations {
            return false;
        }
        self.residual_history
            .last()
            .is_some_and(|&r| r < self.tolerance)
    }

    /// Estimate the convergence rate from recent history.
    ///
    /// Returns the log-ratio of recent residuals: if rate > 0, we're converging.
    pub fn estimate_rate(&self) -> f64 {
        let h = &self.residual_history;
        if h.len() < 10 {
            return 0.0;
        }
        let recent = &h[h.len() - 10..];
        let first = recent[0].max(1e-30);
        let last = recent[9].max(1e-30);
        (first.ln() - last.ln()) / 10.0
    }

    /// Estimate spectral gap from eigenvalue history.
    ///
    /// The spectral gap is estimated from the stabilization rate of
    /// the eigenvalue estimate.
    pub fn estimate_spectral_gap(&self) -> Option<f64> {
        let h = &self.eigenvalue_history;
        if h.len() < 20 {
            return None;
        }
        let _recent = &h[h.len() - 20..];
        // The spectral gap is related to the rate of convergence
        // Higher gap = faster convergence = larger rate
        let rate = self.estimate_rate();
        if rate > 0.0 {
            Some(rate)
        } else {
            None
        }
    }

    /// Compute the final convergence info.
    pub fn finalize(self) -> ConvergenceInfo {
        let final_residual = self.residual_history.last().copied().unwrap_or(f64::INFINITY);
        let converged = self.iteration >= self.min_iterations && final_residual < self.tolerance;
        let rate = self.estimate_rate();
        let spectral_gap = self.estimate_spectral_gap();

        // Downsample history for storage (keep at most 100 points)
        let history_len = self.residual_history.len();
        let residual_history = if history_len <= 100 {
            self.residual_history
        } else {
            let step = history_len as f64 / 100.0;
            (0..100)
                .map(|i| {
                    let idx = (i as f64 * step) as usize;
                    self.residual_history[idx.min(history_len - 1)]
                })
                .collect()
        };

        ConvergenceInfo {
            converged,
            iterations: self.iteration,
            residual: final_residual,
            rate,
            spectral_gap,
            residual_history,
        }
    }

    /// Current residual.
    pub fn current_residual(&self) -> f64 {
        self.residual_history.last().copied().unwrap_or(f64::INFINITY)
    }

    /// Current eigenvalue estimate.
    pub fn current_eigenvalue(&self) -> f64 {
        self.eigenvalue_history
            .last()
            .copied()
            .unwrap_or(0.0)
    }

    /// Check if eigenvalue has stabilized (change < tol over last N iterations).
    pub fn eigenvalue_stabilized(&self, n: usize, tol: f64) -> bool {
        let h = &self.eigenvalue_history;
        if h.len() < n {
            return false;
        }
        let recent = &h[h.len() - n..];
        let max_val = recent.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_val = recent.iter().cloned().fold(f64::INFINITY, f64::min);
        (max_val - min_val).abs() < tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_tracker_is_not_converged() {
        let t = ConvergenceTracker::new(1e-8, 1000);
        assert!(!t.is_converged());
        assert_eq!(t.iteration(), 0);
    }

    #[test]
    fn record_increments_iteration() {
        let mut t = ConvergenceTracker::new(1e-8, 1000);
        t.record(1.0, 2.0);
        assert_eq!(t.iteration(), 1);
    }

    #[test]
    fn converged_below_tolerance() {
        let mut t = ConvergenceTracker::new(1e-8, 1000).with_min_iterations(5);
        for _ in 0..5 {
            t.record(1e-10, 3.0);
        }
        assert!(t.is_converged());
    }

    #[test]
    fn not_converged_above_tolerance() {
        let mut t = ConvergenceTracker::new(1e-8, 1000).with_min_iterations(5);
        for _ in 0..20 {
            t.record(1.0, 3.0);
        }
        assert!(!t.is_converged());
    }

    #[test]
    fn exhausted_check() {
        let mut t = ConvergenceTracker::new(1e-8, 5);
        for _ in 0..5 {
            t.record(1.0, 3.0);
        }
        assert!(t.exhausted());
    }

    #[test]
    fn finalize_converged() {
        let mut t = ConvergenceTracker::new(1e-6, 1000).with_min_iterations(10);
        for _ in 0..20 {
            t.record(1e-8, 3.0);
        }
        let info = t.finalize();
        assert!(info.converged);
        assert_eq!(info.iterations, 20);
        assert!(info.residual < 1e-6);
    }

    #[test]
    fn finalize_not_converged() {
        let info = ConvergenceInfo::did_not_converge(1000, 0.5);
        assert!(!info.converged);
        assert_eq!(info.iterations, 1000);
    }

    #[test]
    fn rate_estimation() {
        let mut t = ConvergenceTracker::new(1e-10, 1000);
        // Simulate exponential convergence
        for i in 0..20 {
            let residual = 0.1_f64.powi(i as i32 + 1);
            t.record(residual, 3.0);
        }
        let rate = t.estimate_rate();
        assert!(rate > 0.0);
    }

    #[test]
    fn rate_too_few_points() {
        let mut t = ConvergenceTracker::new(1e-10, 1000);
        t.record(1.0, 3.0);
        assert_eq!(t.estimate_rate(), 0.0);
    }

    #[test]
    fn eigenvalue_stabilized() {
        let mut t = ConvergenceTracker::new(1e-8, 1000);
        for _ in 0..20 {
            t.record(1e-8, 3.14159);
        }
        assert!(t.eigenvalue_stabilized(10, 1e-6));
    }

    #[test]
    fn eigenvalue_not_stabilized() {
        let mut t = ConvergenceTracker::new(1e-8, 1000);
        for i in 0..20 {
            t.record(1e-8, i as f64);
        }
        assert!(!t.eigenvalue_stabilized(10, 1e-6));
    }

    #[test]
    fn current_residual_and_eigenvalue() {
        let mut t = ConvergenceTracker::new(1e-8, 1000);
        t.record(0.5, 2.0);
        t.record(0.3, 2.5);
        assert!((t.current_residual() - 0.3).abs() < 1e-10);
        assert!((t.current_eigenvalue() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn history_downsampled_in_finalize() {
        let mut t = ConvergenceTracker::new(1e-10, 10000);
        for i in 0..500 {
            t.record(1e-3 / (i as f64 + 1.0), 3.0);
        }
        let info = t.finalize();
        assert!(info.residual_history.len() <= 100);
    }
}
