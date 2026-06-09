//! Tutorial: analog-spectral — eigenvalue estimation as mechanical settling.
//!
//! Run with: `cargo run --example tutorial`

use analog_spectral::{
    AnalogDial, DialBank, SpectralGapAnalysis, SpectralThermostat, Action,
    PrecisionAnalysis,
};

fn main() {
    println!("=== analog-spectral Tutorial ===\n");

    // ---- 1. Single Dial Dynamics ----
    println!("--- 1. Single Dial Dynamics ---");
    let mut dial = AnalogDial::new(3.14, 5.0, 0.5);
    println!("Deadband (spectral gap): {:.4}", dial.deadband());

    dial.position = 0.0;
    dial.velocity = 0.0;
    for _ in 0..500 {
        dial.step(0.01);
    }
    println!("After 500 steps: position = {:.4} (target: 3.14)", dial.position);
    println!("Is settled: {}", dial.is_settled());
    println!();

    // ---- 2. Deadband IS the Spectral Gap ----
    println!("--- 2. Deadband = Spectral Gap ---");
    let tight = AnalogDial::new(5.0, 1.0, 1.0);
    let clear = AnalogDial::new(5.0, 100.0, 0.1);
    println!(
        "Tight gap (low gravity, high friction): deadband = {:.4}",
        tight.deadband()
    );
    println!(
        "Clear gap (high gravity, low friction): deadband = {:.6}",
        clear.deadband()
    );
    println!("Precision: tight = {:.4}, clear = {:.6}", tight.precision(), clear.precision());
    println!();

    // ---- 3. Settle to Convergence ----
    println!("--- 3. Settling to Convergence ---");
    let mut dial = AnalogDial::new(7.5, 20.0, 0.8);
    dial.position = 0.0;
    let steps = dial.settle(0.01, 1e-8);
    println!(
        "Settled in {} steps to {:.8} (target: 7.5)",
        steps, dial.position
    );
    println!();

    // ---- 4. Coupled Dial Bank (Eigenvector Estimation) ----
    println!("--- 4. Coupled Dial Bank ---");
    let coupling = vec![
        vec![2.0, 1.0, 0.0],
        vec![1.0, 3.0, 1.0],
        vec![0.0, 1.0, 2.0],
    ];

    let mut bank = DialBank::new(3, coupling.clone());
    println!("Initial positions: {:?}", bank.read_positions());

    let steps = bank.settle(0.01, 1e-8);
    let positions = bank.read_positions();
    let eigenvalue = bank.eigenvalue_estimate();
    let residual = bank.verify_eigenvector(&coupling);

    println!("Settled in {} steps", steps);
    println!("Positions (eigenvector): [{:.4}, {:.4}, {:.4}]", positions[0], positions[1], positions[2]);
    println!("Eigenvalue (Rayleigh quotient): {:.6}", eigenvalue);
    println!("Residual ||Ax - λx|| / ||x||: {:.2e}", residual);
    println!();

    // ---- 5. Spectral Gap Analysis ----
    println!("--- 5. Spectral Gap Analysis ---");
    let analysis = SpectralGapAnalysis::from_eigenvalues(vec![0.1, 0.5, 1.2, 3.0, 7.5]);
    let (idx, gap) = analysis.largest_gap();
    println!("Largest gap: index {}, value {:.4}", idx, gap);
    println!("Conditioning ratio: {:.4}", analysis.conditioning());
    println!();
    println!("Gap details:");
    for i in 0..analysis.gaps.len() {
        println!(
            "  Gap {}: λ[{}]→λ[{}] = {:.4}, deadband={:.4}, settle_time={:.4}",
            i, i, i + 1, analysis.gaps[i], analysis.deadbands[i], analysis.settling_time(i)
        );
    }
    println!();

    // ---- 6. Spectral Thermostat ----
    println!("--- 6. Spectral Thermostat ---");
    let mut thermo = SpectralThermostat::new(0.5, 0.05);
    println!("Target CR = 0.5, deadband = ±0.05");

    let measurements = [0.42, 0.48, 0.51, 0.55, 0.50, 0.49, 0.50, 0.47];
    for &cr in &measurements {
        let action = thermo.measure(cr);
        let action_str = match action {
            Action::IncreaseCR => "HEAT ↑",
            Action::DecreaseCR => "COOL ↓",
            Action::DoNothing => "HOLD ─",
        };
        println!("  CR={:.2} → {} ({:?})", cr, action_str, thermo.state());
    }
    println!("Hysteresis: {:.2}", thermo.hysteresis());
    println!("Stable for {} consecutive measurements", thermo.stability_duration());
    println!();

    // ---- 7. Precision Analysis ----
    println!("--- 7. Precision Analysis ---");
    let dial = AnalogDial::new(10.0, 50.0, 0.1);
    let analysis = PrecisionAnalysis::for_dial(&dial);
    println!("Thermal noise: {:.2e}", analysis.thermal_noise);
    println!("Friction noise: {:.2e}", analysis.friction_noise);
    println!("Gravity precision: {:.2e}", analysis.gravity_precision);
    println!("Effective bit depth: {:.1} (vs 53 for f64)", analysis.eigenvalue_bits());
    println!("Precision ratio: {:.4}", analysis.vs_digital());
    println!("Dominant noise: {}", analysis.dominant_noise_source());

    // Eigenvector precision with coupling matrix
    let coupling = vec![
        vec![2.0, 1.0, 0.0],
        vec![1.0, 3.0, 1.0],
        vec![0.0, 1.0, 2.0],
    ];
    println!(
        "Eigenvector bits (3×3 matrix): {:.1}",
        analysis.eigenvector_bits(&coupling)
    );

    println!("\n=== Tutorial Complete ===");
}
