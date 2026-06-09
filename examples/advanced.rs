//! Advanced: Full spectral analysis pipeline — from graph to eigenvalues via analog dials.
//!
//! Demonstrates computing eigenvalues of graph Laplacians using analog dial banks,
//! analyzing spectral gaps, and controlling computation with a spectral thermostat.
//!
//! Run with: `cargo run --example advanced`

use analog_spectral::{
    AnalogDial, DialBank, SpectralGapAnalysis, SpectralThermostat, Action,
    PrecisionAnalysis,
};

fn main() {
    println!("=== Advanced: Full Spectral Analysis Pipeline ===\n");

    // ---- Build several graph Laplacians and analyze spectra ----

    // Path graph Laplacian (5 nodes)
    let path_adj = vec![
        vec![0.0, 1.0, 0.0, 0.0, 0.0],
        vec![1.0, 0.0, 1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0, 0.0, 1.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0],
    ];

    // Compute Laplacian L = D - A
    let laplacian = compute_laplacian(&path_adj);

    println!("Path Graph Laplacian (5 nodes):");
    print_matrix(&laplacian);
    println!();

    // ---- Estimate dominant eigenvalue via dial bank ----
    println!("--- Eigenvalue Estimation via Dial Bank ---");
    let mut bank = DialBank::new(5, laplacian.clone());
    let steps = bank.settle(0.01, 1e-8);
    let positions = bank.read_positions();
    let eigenvalue = bank.eigenvalue_estimate();
    let residual = bank.verify_eigenvector(&laplacian);

    println!("Settled in {} steps", steps);
    println!("Eigenvector estimate: [{:.4}, {:.4}, {:.4}, {:.4}, {:.4}]",
        positions[0], positions[1], positions[2], positions[3], positions[4]);
    println!("Eigenvalue (Rayleigh quotient): {:.6}", eigenvalue);
    println!("Residual: {:.2e}", residual);
    println!();

    // ---- Spectral gap analysis ----
    // Use known path graph eigenvalues for analysis
    // For path graph of size 5: eigenvalues ≈ 0, 0.382, 1.0, 2.0, 3.618
    let path_eigenvalues = vec![0.0, 0.382, 1.0, 2.0, 3.618];
    println!("--- Spectral Gap Analysis ---");
    let gap_analysis = SpectralGapAnalysis::from_eigenvalues(path_eigenvalues.clone());
    let (max_idx, max_gap) = gap_analysis.largest_gap();
    println!("Eigenvalues: {:?}", path_eigenvalues.iter().map(|v| format!("{:.3}", v)).collect::<Vec<_>>());
    println!("Largest gap: index {} (between λ[{}] and λ[{}]), width = {:.4}",
        max_idx, max_idx, max_idx + 1, max_gap);
    println!("Conditioning ratio: {:.4}", gap_analysis.conditioning());
    println!();

    println!("Gap details:");
    for i in 0..gap_analysis.gaps.len() {
        println!(
            "  λ[{}]→λ[{}]: gap={:.4}, deadband={:.4}, settle_time={:.4}",
            i, i + 1, gap_analysis.gaps[i], gap_analysis.deadbands[i],
            gap_analysis.settling_time(i)
        );
    }
    println!();

    // ---- Compare different gravity/friction regimes ----
    println!("--- Gravity/Friction Regimes ---");
    let regimes = [
        ("High precision", 100.0, 0.01),
        ("Standard", 10.0, 0.5),
        ("Low precision", 1.0, 2.0),
        ("Ultra-tight", 1000.0, 0.001),
    ];

    for (name, gravity, friction) in &regimes {
        let mut dial = AnalogDial::new(5.0, *gravity, *friction);
        let deadband = dial.deadband();
        dial.position = 0.0;
        let steps = dial.settle(0.01, 1e-10);
        let precision = PrecisionAnalysis::for_dial(&dial);
        println!(
            "  {}: gravity={}, friction={}, deadband={:.6}, steps={}, bits={:.1}, noise={}",
            name, gravity, friction, deadband, steps,
            precision.eigenvalue_bits(),
            precision.dominant_noise_source()
        );
    }
    println!();

    // ---- Spectral thermostat: controlling conservation ratio ----
    println!("--- Thermostat Control Sequence ---");
    let mut thermo = SpectralThermostat::new(0.5, 0.05);

    // Simulate a control loop: CR fluctuates, thermostat decides
    let cr_readings = [0.30, 0.35, 0.42, 0.48, 0.49, 0.50, 0.51, 0.53, 0.56, 0.52, 0.50, 0.49];
    for cr in &cr_readings {
        let action = thermo.measure(*cr);
        let indicator = match action {
            Action::IncreaseCR => "↑ HEAT",
            Action::DecreaseCR => "↓ COOL",
            Action::DoNothing => "● HOLD",
        };
        println!("  CR={:.2} {}", cr, indicator);
    }
    println!("Hysteresis ratio: {:.2} (fraction of measurements that triggered action)", thermo.hysteresis());
    println!("Stable for {} consecutive readings", thermo.stability_duration());
    println!();

    // ---- Precision comparison across dial configurations ----
    println!("--- Precision vs Digital f64 ---");
    let configs = [
        ("Settle=1.0, g=10, f=0.5", 1.0, 10.0, 0.5),
        ("Settle=10.0, g=10, f=0.5", 10.0, 10.0, 0.5),
        ("Settle=100.0, g=50, f=0.1", 100.0, 50.0, 0.1),
    ];

    for (name, setpoint, gravity, friction) in &configs {
        let dial = AnalogDial::new(*setpoint, *gravity, *friction);
        let analysis = PrecisionAnalysis::for_dial(&dial);
        println!(
            "  {}: {:.1} bits ({:.2} of f64), dominant={}",
            name, analysis.eigenvalue_bits(), analysis.vs_digital(),
            analysis.dominant_noise_source()
        );
    }

    println!("\n=== Analysis Complete ===");
}

fn compute_laplacian(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = adj.len();
    let mut lap = vec![vec![0.0; n]; n];
    for i in 0..n {
        let degree: f64 = adj[i].iter().sum();
        lap[i][i] = degree;
        for j in 0..n {
            if i != j {
                lap[i][j] = -adj[i][j];
            }
        }
    }
    lap
}

fn print_matrix(m: &[Vec<f64>]) {
    for row in m {
        let s: Vec<String> = row.iter().map(|v| format!("{:6.2}", v)).collect();
        println!("  [{}]", s.join(", "));
    }
}
