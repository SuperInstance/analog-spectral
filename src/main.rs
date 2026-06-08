//! Analog spectral: physical dials, gravity deadbands, spectral computation.

use analog_spectral::{AnalogDial, SpectralGapAnalysis};

fn main() {
    println!("analog-spectral: physical dials, gravity deadbands, spectral computation");

    // Demo: single dial settling
    let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
    dial.position = 0.0;
    let steps = dial.settle(0.01, 1e-10);
    println!("Dial settled in {} steps to {:.4} (setpoint=5.0)", steps, dial.position);

    // Demo: spectral gap analysis
    let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);
    let (idx, gap) = analysis.largest_gap();
    println!("Largest spectral gap: index={}, value={:.4}", idx, gap);
}
