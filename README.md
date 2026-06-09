# analog-spectral

**Eigenvalue estimation as mechanical settling — physical dials converge under gravity, friction creates deadbands equal to spectral gaps, and the thermostat IS the algorithm.**

> *What if eigenvalue computation isn't a numerical procedure but a physical process? Dials settling under gravity toward eigenvalue setpoints, with friction creating the deadbands that are spectral gaps. The Jacobi algorithm is just a bunch of dials finding equilibrium.*

---

## The Problem

Eigenvalue computation is fundamental to spectral graph theory, conservation analysis, and the entire SuperInstance ecosystem. But numerical algorithms abstract away the physics. What if you could *see* the computation happening — dials turning, springs pulling, friction resisting — and the settled positions *were* the eigenvalues?

This matters because the physical analogy reveals structure that pure numerics hides: the deadband (friction/gravity) is exactly the spectral gap between eigenvalues, the settling time predicts convergence rate, and the precision is bounded by thermal noise in exactly the way floating-point is bounded by mantissa bits.

## The Key Insight

**A damped harmonic oscillator converging to a setpoint IS eigenvalue iteration.**

Consider a single analog dial:
- **Position** = current eigenvalue estimate
- **Setpoint** = true eigenvalue
- **Gravity** = restoring force strength (convergence rate)
- **Friction** = damping coefficient
- **Deadband** = friction/gravity = spectral gap

The equation of motion is:
```
m·ẍ = −gravity·(x − setpoint) − friction·ẋ
```

Within the deadband (|x − setpoint| < friction/gravity), friction dominates and the dial stops. This deadband IS the spectral gap — the region where the eigenvalue estimate is "close enough" that numerical noise dominates further convergence.

For N coupled dials, the coupling matrix determines inter-dial forces. When the bank settles, the dial positions approximate eigenvectors and the Rayleigh quotient gives the eigenvalue. **This IS the Jacobi eigenvalue algorithm, but phrased as physical dynamics.**

## Architecture

```
 ┌─────────────────────────────────────────────────┐
 │                    DialBank                      │
 │  N coupled dials + coupling matrix              │
 │                                                 │
 │  ┌─────────┐ ┌─────────┐       ┌─────────┐    │
 │  │ Dial 0  │ │ Dial 1  │ ···   │ Dial N  │    │
 │  │ pos=λ₀  │ │ pos=λ₁  │       │ pos=λₙ  │    │
 │  └────┬────┘ └────┬────┘       └────┬────┘    │
 │       │            │                 │          │
 │       └────────────┼─────────────────┘          │
 │                    │                            │
 │            coupling matrix A                    │
 │         (inter-dial forces)                     │
 └────────────┬────────────────────────────────────┘
              │
              ▼
 ┌─────────────────────────┐    ┌──────────────────────┐
 │    SpectralGapAnalysis  │    │  PrecisionAnalysis   │
 │                         │    │                      │
 │  eigenvalues → gaps     │    │  thermal noise       │
 │  gaps → deadbands       │    │  friction noise      │
 │  gaps → settling times  │    │  effective bit depth │
 └─────────────────────────┘    │  vs digital f64      │
                                └──────────────────────┘
              │
              ▼
 ┌─────────────────────────┐
 │  SpectralThermostat     │
 │                         │
 │  measure(CR) → Action   │
 │  ┌─────────────────┐   │
 │  │ within deadband │──▶─┼──▶ DoNothing (stable)
 │  │ below setpoint  │──▶─┼──▶ IncreaseCR (heating)
 │  │ above setpoint  │──▶─┼──▶ DecreaseCR (cooling)
 │  └─────────────────┘   │
 │  + hysteresis tracking  │
 └─────────────────────────┘
```

## Quick Start

```rust
use analog_spectral::{AnalogDial, SpectralThermostat, Action};

// A dial seeking eigenvalue ≈ 5.0
// gravity=10, friction=1 → deadband = 0.1
let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
assert!((dial.deadband() - 0.1).abs() < 1e-12);

// Displace it and let it settle
dial.position = 0.0;
let steps = dial.settle(0.01, 1e-6);
println!("Settled in {steps} steps to {:.6}", dial.position);

// Spectral thermostat for conservation ratio control
let mut thermo = SpectralThermostat::new(0.5, 0.05);
assert_eq!(thermo.measure(0.48), Action::IncreaseCR); // below setpoint
assert_eq!(thermo.measure(0.55), Action::DecreaseCR); // above setpoint
assert_eq!(thermo.measure(0.52), Action::DoNothing);   // within deadband
```

## Tutorial

### 1. Single Dial Dynamics

A dial under gravity and friction converges to its setpoint:

```rust
use analog_spectral::AnalogDial;

let mut dial = AnalogDial::new(3.14, 5.0, 0.5);
dial.position = 0.0; // start far from setpoint
dial.velocity = 0.0;

// Step it forward
for _ in 0..1000 {
    dial.step(0.01);
}
println!("Position: {:.4} (target: 3.14)", dial.position);

// Or use settle() to iterate until convergence
dial.position = 0.0;
let steps = dial.settle(0.01, 1e-10);
println!("Settled in {steps} steps");
```

### 2. The Deadband IS the Spectral Gap

The deadband width = friction/gravity. This is the physical analog of the spectral gap:

```rust
use analog_spectral::AnalogDial;

// Tight spectral gap: low gravity, high friction → wide deadband → slow settling
let mut tight = AnalogDial::new(5.0, 1.0, 1.0);
assert!((tight.deadband() - 1.0).abs() < 1e-12); // deadband = 1.0

// Clear spectral gap: high gravity, low friction → narrow deadband → fast settling
let mut clear = AnalogDial::new(5.0, 100.0, 0.1);
assert!((clear.deadband() - 0.001).abs() < 1e-12); // deadband = 0.001
```

### 3. Coupled Dial Banks (Eigenvector Estimation)

N dials coupled by a matrix approximate eigenvectors:

```rust
use analog_spectral::DialBank;

// A 3×3 coupling matrix
let coupling = vec![
    vec![2.0, 1.0, 0.0],
    vec![1.0, 3.0, 1.0],
    vec![0.0, 1.0, 2.0],
];

let mut bank = DialBank::new(3, coupling.clone());
let steps = bank.settle(0.01, 1e-8);
println!("Bank settled in {steps} steps");

// Read eigenvector components and eigenvalue estimate
let positions = bank.read_positions();
let eigenvalue = bank.eigenvalue_estimate();
println!("Eigenvector estimate: {:?}", positions);
println!("Eigenvalue estimate (Rayleigh quotient): {:.6}", eigenvalue);

// Verify quality via residual ||Ax − λx|| / ||x||
let residual = bank.verify_eigenvector(&coupling);
println!("Residual: {:.2e}", residual);
```

### 4. Spectral Gap Analysis

Given eigenvalues, compute gaps, deadbands, and settling times:

```rust
use analog_spectral::SpectralGapAnalysis;

let analysis = SpectralGapAnalysis::from_eigenvalues(vec![0.1, 0.5, 1.2, 3.0, 7.5]);

let (idx, gap) = analysis.largest_gap();
println!("Largest gap: index {}, value {:.4}", idx, gap);

println!("Conditioning ratio: {:.4}", analysis.conditioning());
for i in 0..analysis.gaps.len() {
    println!(
        "Gap {}: width={:.4}, deadband={:.4}, settling_time={:.4}",
        i, analysis.gaps[i], analysis.deadbands[i], analysis.settling_time(i)
    );
}
```

### 5. Spectral Thermostat

Deadband-based control for spectral quantities:

```rust
use analog_spectral::{SpectralThermostat, Action};

let mut thermo = SpectralThermostat::new(0.5, 0.05); // target CR=0.5, deadband=0.05

// Simulate fluctuating measurements
let measurements = [0.42, 0.48, 0.51, 0.55, 0.50, 0.49, 0.50];
for cr in &measurements {
    let action = thermo.measure(*cr);
    println!("CR={:.2} → {:?}", cr, action);
}

println!("Hysteresis: {:.2}", thermo.hysteresis()); // fraction of actions
println!("Stable for {} consecutive measurements", thermo.stability_duration());
```

### 6. Precision Analysis

Compare analog precision to digital f64:

```rust
use analog_spectral::{AnalogDial, PrecisionAnalysis};

let dial = AnalogDial::new(10.0, 50.0, 0.1);
let analysis = PrecisionAnalysis::for_dial(&dial);

println!("Effective bits: {:.1} (vs 53 for f64)", analysis.eigenvalue_bits());
println!("Precision ratio: {:.4}", analysis.vs_digital());
println!("Dominant noise: {}", analysis.dominant_noise_source());
```

## API Reference

| Type | Module | Key Methods |
|------|--------|-------------|
| `AnalogDial` | `analog_dial` | `new()`, `step()`, `settle()`, `deadband()`, `is_settled()`, `precision()` |
| `DialBank` | `dial_bank` | `new()`, `step()`, `settle()`, `read_positions()`, `eigenvalue_estimate()`, `verify_eigenvector()` |
| `SpectralGapAnalysis` | `spectral_gap` | `from_eigenvalues()`, `largest_gap()`, `deadband_at()`, `conditioning()`, `settling_time()` |
| `SpectralThermostat` | `thermostat` | `new()`, `measure()`, `stability_duration()`, `hysteresis()`, `state()` |
| `PrecisionAnalysis` | `precision` | `for_dial()`, `eigenvalue_bits()`, `eigenvector_bits()`, `dominant_noise_source()`, `vs_digital()` |
| `ThermostatState` | `thermostat` | `Heating`, `Cooling`, `Stable` |
| `Action` | `thermostat` | `IncreaseCR`, `DecreaseCR`, `DoNothing` |

## Ecosystem Role

`analog-spectral` is the **physical analog layer** of the SuperInstance ecosystem:

- **[constraint-theory-core](https://github.com/SuperInstance/constraint-theory-core)** — Uses deadbands and spectral gaps for constraint analysis
- **[dial-ecology](https://github.com/SuperInstance/dial-ecology)** — Lotka-Volterra competition dynamics for traditions
- **[wave-conservation](https://github.com/SuperInstance/wave-conservation)** — Wave propagation and spectral analysis on graphs
- **analog-spectral** — Physical dials, deadbands, and eigenvalue estimation (this crate)

The deadband concept (friction/gravity = spectral gap) unifies all crates. The thermostat controller provides the decision logic for conservation ratio management across the ecosystem.

## Installation

```toml
[dependencies]
analog-spectral = "0.1"
```

## License

MIT
