use crate::analog_dial::AnalogDial;

/// How precise can analog computation be?
/// Analog precision limited by: thermal noise, friction, gravity fluctuations.
pub struct PrecisionAnalysis {
    pub thermal_noise: f64,
    pub friction_noise: f64,
    pub gravity_precision: f64,
    pub effective_bits: f64,
}

impl PrecisionAnalysis {
    /// Analyze precision of analog dial computation.
    pub fn for_dial(dial: &AnalogDial) -> PrecisionAnalysis {
        // Thermal noise: kT at room temperature, scaled to dial range
        let k_boltzmann = 1.38064852e-23;
        let temperature = 300.0; // room temp in Kelvin
        let thermal = k_boltzmann * temperature;

        // Friction noise: stick-slip proportional to friction coefficient
        let friction_noise = dial.friction * 1e-3;

        // Gravity precision: how well we know gravity (~1e-9 relative)
        let gravity_precision = dial.gravity * 1e-9;

        // Total noise floor
        let noise_floor = thermal + friction_noise + gravity_precision;

        // Dynamic range: setpoint span
        let dynamic_range = dial.setpoint.abs().max(dial.deadband());

        // Effective bits = log2(dynamic_range / noise_floor)
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

    /// Bits of eigenvalue precision.
    pub fn eigenvalue_bits(&self) -> f64 {
        self.effective_bits
    }

    /// Bits of eigenvector precision (degraded by coupling condition).
    pub fn eigenvector_bits(&self, coupling: &[Vec<f64>]) -> f64 {
        // Condition number estimate from coupling matrix
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

        // Rough: eigenvector bits = eigenvalue bits - log2(condition_estimate)
        let cond_estimate = max_entry * (n as f64);
        let cond_loss = cond_estimate.log2().max(0.0);
        (self.effective_bits - cond_loss).max(1.0)
    }

    /// What degrades first?
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
    pub fn vs_digital(&self) -> f64 {
        self.effective_bits / 53.0
    }
}
