//! Precision analysis for analog eigenvalue computation.
//!
//! Analog precision is limited by thermal noise, friction stick-slip,
//! and gravity fluctuations. This module quantifies the effective bit
//! depth of analog dial computation and compares it to digital f64.

use crate::analog_dial::AnalogDial;

/// Precision analysis comparing analog computation to digital limits.
///
/// Estimates the effective bit depth of eigenvalue and eigenvector
/// estimates produced by analog dials, accounting for physical noise sources.
pub struct PrecisionAnalysis {
    /// Thermal noise floor (kT at room temperature, scaled).
    pub thermal_noise: f64,
    /// Friction-induced stick-slip noise.
    pub friction_noise: f64,
    /// Gravity fluctuation precision limit.
    pub gravity_precision: f64,
    /// Effective bit depth of the eigenvalue estimate.
    pub effective_bits: f64,
}

impl PrecisionAnalysis {
    /// Analyze the precision limits of an analog dial.
    ///
    /// Computes noise from three sources:
    /// - **Thermal**: kT at 300K, scaled to dial range
    /// - **Friction**: stick-slip proportional to friction coefficient
    /// - **Gravity**: drift from gravity uncertainty (~10⁻⁹ relative)
    pub fn for_dial(dial: &AnalogDial) -> PrecisionAnalysis {
        let k_boltzmann = 1.38064852e-23;
        let temperature = 300.0;
        let thermal = k_boltzmann * temperature;

        let friction_noise = dial.friction * 1e-3;
        let gravity_precision = dial.gravity * 1e-9;

        let noise_floor = thermal + friction_noise + gravity_precision;
        let dynamic_range = dial.setpoint.abs().max(dial.deadband());

        let effective_bits = if noise_floor > 0.0 && dynamic_range > 0.0 {
            (dynamic_range / noise_floor).log2()
        } else {
            0.0
        };

        PrecisionAnalysis {
            thermal_noise: thermal,
            friction_noise,
            gravity_precision,
            effective_bits: effective_bits.max(0.0),
        }
    }

    /// Effective bits of eigenvalue precision.
    ///
    /// Analog dials achieve roughly log₂(dynamic_range / noise_floor) bits.
    pub fn eigenvalue_bits(&self) -> f64 {
        self.effective_bits
    }

    /// Effective bits of eigenvector precision, degraded by coupling condition.
    ///
    /// Eigenvectors are more sensitive to noise than eigenvalues when
    /// the coupling matrix is ill-conditioned.
    pub fn eigenvector_bits(&self, coupling: &[Vec<f64>]) -> f64 {
        let n = coupling.len();
        if n == 0 {
            return 0.0;
        }

        let mut max_entry = 0.0f64;
        for row in coupling {
            for &v in row {
                max_entry = max_entry.max(v.abs());
            }
        }
        if max_entry == 0.0 {
            return 0.0;
        }

        let cond_estimate = max_entry * (n as f64);
        let cond_loss = cond_estimate.log2().max(0.0);
        (self.effective_bits - cond_loss).max(1.0)
    }

    /// Identify the dominant noise source limiting precision.
    ///
    /// Returns one of: "friction (stick-slip)", "gravity drift", or "thermal noise (kT)".
    pub fn dominant_noise_source(&self) -> String {
        if self.friction_noise >= self.thermal_noise && self.friction_noise >= self.gravity_precision {
            "friction (stick-slip)".to_string()
        } else if self.gravity_precision >= self.thermal_noise {
            "gravity drift".to_string()
        } else {
            "thermal noise (kT)".to_string()
        }
    }

    /// Compare analog precision to digital f64 (53-bit mantissa).
    ///
    /// Returns a ratio: 1.0 means analog matches f64 precision,
    /// values below 1.0 indicate how much precision is lost.
    pub fn vs_digital(&self) -> f64 {
        self.effective_bits / 53.0
    }
}
