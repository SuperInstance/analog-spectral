# analog-spectral

Physical analog dials and spectral thermostats — eigenvalue estimation as mechanical settling, deadbands as spectral gaps, and thermostat logic for conservation control. Rust, 25 tests.

## What This Gives You

- **Analog dials** — physical simulation of a dial under gravity settling to an eigenvalue estimate
- **Deadband = spectral gap** — the friction/gravity ratio gives the settling deadband
- **Dial banks** — arrays of dials estimating an entire eigenvalue spectrum
- **Spectral thermostats** — deadband logic to control conservation ratios (increase/decrease/hold)
- **Precision analysis** — relate deadband width to measurement precision
- **25 tests** — verified physics, convergence, and thermostat behavior

## Quick Start

```rust
use analog_spectral::{AnalogDial, SpectralThermostat, Action};

// A dial settling to eigenvalue ≈ 5.0
let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
// gravity=10, friction=1 → deadband = 0.1
assert!((dial.deadband() - 0.1).abs() < 1e-12);

// Settle from a displaced position
dial.position = 0.0;
let steps = dial.settle(0.01, 1e-6);
println!("Settled in {steps} steps");

// Spectral thermostat — control conservation ratio
let mut thermo = SpectralThermostat::new(0.5, 0.05); // target CR=0.5, deadband=0.05
let action = thermo.measure(0.48); // current CR
assert_eq!(action, Action::IncreaseCR); // below setpoint
```

## The Physics

An analog dial has:
- **Setpoint** = target eigenvalue
- **Gravity** = restoring force strength
- **Friction** = damping
- **Deadband** = friction/gravity = minimum resolvable distance (the spectral gap)
- **Settling time** ∝ 1/gravity (higher gravity → faster convergence)

A spectral thermostat applies deadband logic to conservation ratio control:
- Within deadband → DoNothing (stable)
- Below setpoint → IncreaseCR (heating)
- Above setpoint → DecreaseCR (cooling)

## API Reference

| Type | Description |
|---|---|
| `AnalogDial` | Physical dial: setpoint, gravity, friction, deadband |
| `DialBank` | Array of dials for full spectrum estimation |
| `SpectralThermostat` | Deadband controller for conservation ratio |
| `SpectralGapAnalysis` | Eigenvalue gap analysis via dial settling |
| `PrecisionAnalysis` | Deadband-to-precision relationship |

## How It Fits

The **analog physics layer** of the constraint theory ecosystem:

- [analog-spline-theory](https://github.com/SuperInstance/analog-spline-theory) — theoretical foundations
- [constraint-theory-core](https://github.com/SuperInstance/constraint-theory-core) — uses deadbands and spectral gaps
- [conservation-spectral-python](https://github.com/SuperInstance/conservation-spectral-python) — spectral analysis in Python

## Testing

```bash
cargo test  # 25 tests
```

## Installation

```bash
cargo add analog-spectral
```

## License

MIT
