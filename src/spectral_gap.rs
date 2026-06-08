//! Spectral gap analysis — eigenvalue gaps as deadband widths.
//!
//! The spectral gap between consecutive eigenvalues determines the deadband
//! of the analog system. A large gap means a wide deadband and fast settling
//! (clear mode separation). A small gap means narrow deadband and slow
//! settling (nearly degenerate modes).

/// Analysis of eigenvalue gaps and their physical deadband equivalents.
///
/// Given a sorted list of eigenvalues, computes the gaps between consecutive
/// values and normalizes them into deadband widths relative to the largest
/// eigenvalue.
pub struct SpectralGapAnalysis {
    /// Sorted eigenvalues (ascending).
    pub eigenvalues: Vec<f64>,
    /// Gaps between consecutive eigenvalues: gaps[i] = λ[i+1] − λ[i].
    pub gaps: Vec<f64>,
    /// Normalized deadbands: deadbands[i] = gaps[i] / max(λ).
    pub deadbands: Vec<f64>,
}

impl SpectralGapAnalysis {
    /// Compute spectral gap analysis from a set of eigenvalues.
    ///
    /// Eigenvalues are sorted ascending. Gaps and deadbands are computed
    /// for each consecutive pair.
    pub fn from_eigenvalues(mut eigenvalues: Vec<f64>) -> SpectralGapAnalysis {
        eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let max_ev = eigenvalues.iter().cloned().fold(0.0f64, f64::max).max(1e-15);

        let mut gaps = Vec::new();
        let mut deadbands = Vec::new();

        for i in 0..eigenvalues.len().saturating_sub(1) {
            let gap = eigenvalues[i + 1] - eigenvalues[i];
            gaps.push(gap);
            deadbands.push(gap / max_ev);
        }

        SpectralGapAnalysis {
            eigenvalues,
            gaps,
            deadbands,
        }
    }

    /// Find the largest spectral gap — the most stable deadband.
    ///
    /// Returns `(gap_index, gap_value)`. The largest gap indicates where
    /// the spectrum has the clearest separation between mode clusters.
    pub fn largest_gap(&self) -> (usize, f64) {
        if self.gaps.is_empty() {
            return (0, 0.0);
        }
        let mut best_idx = 0;
        let mut best_gap = self.gaps[0];
        for (i, &g) in self.gaps.iter().enumerate() {
            if g > best_gap {
                best_gap = g;
                best_idx = i;
            }
        }
        (best_idx, best_gap)
    }

    /// Deadband width at a specific gap index, normalized by max eigenvalue.
    pub fn deadband_at(&self, index: usize) -> f64 {
        if index < self.deadbands.len() {
            self.deadbands[index]
        } else {
            0.0
        }
    }

    /// Conditioning ratio: min_gap / max_gap.
    ///
    /// Values near 1.0 indicate well-separated eigenvalues.
    /// Values near 0.0 indicate nearly-degenerate modes.
    pub fn conditioning(&self) -> f64 {
        if self.gaps.is_empty() {
            return 0.0;
        }
        let min_gap = self.gaps.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_gap = self.gaps.iter().cloned().fold(0.0f64, f64::max);
        if max_gap > 0.0 {
            min_gap / max_gap
        } else {
            0.0
        }
    }

    /// Estimated settling time for a given gap, proportional to 1/gap.
    ///
    /// Larger gaps settle faster. Returns infinity for zero-width gaps
    /// (degenerate eigenvalues).
    pub fn settling_time(&self, index: usize) -> f64 {
        if index < self.gaps.len() && self.gaps[index] > 0.0 {
            1.0 / self.gaps[index]
        } else {
            f64::INFINITY
        }
    }
}
