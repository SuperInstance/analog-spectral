/// The spectral gap determines the deadband of the system.
/// Large gap → wide deadband → fast settling (clear separation).
/// Small gap → narrow deadband → slow settling (degenerate modes).
pub struct SpectralGapAnalysis {
    pub eigenvalues: Vec<f64>,
    pub gaps: Vec<f64>,
    pub deadbands: Vec<f64>,
}

impl SpectralGapAnalysis {
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

    /// The largest gap = the most stable deadband.
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

    /// Deadband at a specific gap index.
    pub fn deadband_at(&self, index: usize) -> f64 {
        if index < self.deadbands.len() {
            self.deadbands[index]
        } else {
            0.0
        }
    }

    /// Conditioning: min_gap / max_gap.
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

    /// Settling time estimate: proportional to 1/gap.
    pub fn settling_time(&self, index: usize) -> f64 {
        if index < self.gaps.len() && self.gaps[index] > 0.0 {
            1.0 / self.gaps[index]
        } else {
            f64::INFINITY
        }
    }
}
