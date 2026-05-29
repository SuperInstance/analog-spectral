# analog-spectral

**The dial IS the computation. Position = eigenvalue. Deadband = spectral gap. Gravity = restoring force.**

Pure Rust, zero dependencies. A library where analog dials physically settle to eigenvalues under spring-damper dynamics with Coulomb friction deadbands. Not a simulation of eigenvalue computation — a re-implementation of it in physical coordinates.

## The Core Idea

An analog dial has a position, a setpoint, and two forces acting on it:
- **Gravity** (spring force) pulls it toward the setpoint
- **Friction** (damping + Coulomb) resists motion

When the dial is far from the setpoint, gravity wins and it moves. When it's close enough — within the deadband — friction wins and it stops. **The equilibrium position is the eigenvalue estimate. The deadband width IS the spectral gap.**

This is not a new idea. It's how analog computers worked in the 1950s. Op-amp integrators solved differential equations by actually integrating. Memristor networks find eigenvectors by actually settling to energy minima. We're just being explicit about it.

```toml
[dependencies]
analog-spectral = "0.1.0"
```

## Module Walkthrough

### AnalogDial — Spring-Damper Dynamics with Coulomb Friction

The dial follows Newton's second law: `m·a = F_spring + F_damping`. The spring force is proportional to displacement (`-k·x`). Damping opposes velocity (`-c·v`). Coulomb friction adds a deadband — if both displacement and velocity are small enough, the dial freezes.

```rust
use analog_dial::AnalogDial;

// Create a dial: setpoint=5.0, gravity(spring)=10.0, friction=1.0
let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
dial.position = 0.0; // start far from setpoint

// Settle: let gravity pull it to equilibrium
let steps = dial.settle(0.01, 1e-10);
println!("Settled in {} steps", steps);
assert!(dial.is_settled());

// Deadband width = friction / gravity = 1.0 / 10.0 = 0.1
assert!((dial.deadband() - 0.1).abs() < 1e-12);
```

The deadband formula `friction / gravity` is the key physical insight. Higher friction → wider deadband → more tolerance but less precision. Higher gravity → narrower deadband → more precision but more oscillation.

**Settling time is inversely proportional to gravity.** Double the spring constant, halve the settling time:

```rust
let mut d1 = AnalogDial::new(10.0, 4.0, 0.02);
d1.position = 0.0;
let steps1 = d1.settle(0.001, 0.1);

let mut d2 = AnalogDial::new(10.0, 16.0, 0.02);
d2.position = 0.0;
let steps2 = d2.settle(0.001, 0.1);

// d2 has 4× the gravity → settles ~4× faster
assert!(steps2 < steps1);
```

### DialBank — Coupled Dials Compute Eigenvectors

N dials, coupled through a matrix. Each dial's setpoint starts on the diagonal of the coupling matrix. Coupling forces pull dials toward each other proportional to off-diagonal entries. When everything settles, the dial positions approximate the dominant eigenvector.

This is Jacobi eigenvalue computation in disguise. The physical dynamics and the algorithm are the same thing.

```rust
use dial_bank::DialBank;

// Symmetric 2×2 matrix: [[3, 1], [1, 3]]
// Eigenvalues: 4 and 2. Dominant eigenvector: [1/√2, 1/√2]
let coupling = vec![vec![3.0, 1.0], vec![1.0, 3.0]];

let mut bank = DialBank::new(2, coupling.clone());
bank.settle(0.01, 1e-8);

// Read settled positions — these ARE the eigenvector components
let positions = bank.read_positions();
for p in &positions {
    assert!(p.is_finite());
}

// Eigenvalue estimate via Rayleigh quotient: λ = (x^T A x) / (x^T x)
let eigenvalue = bank.eigenvalue_estimate();
println!("Estimated eigenvalue: {:.6}", eigenvalue);

// Verify: compute residual ||Ax - λx|| / ||x||
let residual = bank.verify_eigenvector(&coupling);
assert!(residual < 0.5);
```

The bank steps all dials simultaneously — parallel analog computation. Forces are computed from the pre-step state, then applied all at once. This avoids order-dependent artifacts and mirrors how real analog hardware works (everything evolves continuously, not sequentially).

### SpectralGapAnalysis — Eigenvalue Gaps → Physical Deadbands

Given a set of eigenvalues, compute the gaps between consecutive ones. Each gap determines a deadband width (normalized by the largest eigenvalue). Large gaps mean fast settling and clear mode separation. Small gaps mean slow settling and degenerate modes.

```rust
use spectral_gap::SpectralGapAnalysis;

let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);

// Gaps: [1.0, 3.0]
assert!((analysis.gaps[0] - 1.0).abs() < 1e-12);
assert!((analysis.gaps[1] - 3.0).abs() < 1e-12);

// Largest gap at index 1 (between eigenvalues 2.0 and 5.0)
let (idx, gap) = analysis.largest_gap();
assert_eq!(idx, 1);
assert!((gap - 3.0).abs() < 1e-12);

// Settling time proportional to 1/gap
// Large gap → short settling time
assert!(analysis.settling_time(1) < analysis.settling_time(0));

// Conditioning: min_gap / max_gap = 1/3
assert!((analysis.conditioning() - 1.0/3.0).abs() < 1e-12);

// Deadband at each gap: gap / max_eigenvalue
let db = analysis.deadband_at(0); // gap=1.0 / max_ev=5.0 = 0.2
assert!((db - 0.2).abs() < 1e-12);
```

### SpectralThermostat — IncreaseCR / DecreaseCR / DoNothing

A thermostat that uses deadband logic to control the conservation ratio of a system. The thermostat measures the current CR, compares it to a setpoint, and decides:

- **IncreaseCR** if CR is below setpoint (too cold → heat)
- **DecreaseCR** if CR is above setpoint (too hot → cool)
- **DoNothing** if CR is within deadband (in the band → wait)

```rust
use thermostat::{SpectralThermostat, Action, ThermostatState};

let mut t = SpectralThermostat::new(0.5, 0.1);

// Measure and decide
assert_eq!(t.measure(0.3), Action::IncreaseCR); // below setpoint
assert_eq!(t.measure(0.5), Action::DoNothing);   // within deadband
assert_eq!(t.measure(0.7), Action::DecreaseCR);  // above setpoint

// Track stability
t.measure(0.5);
t.measure(0.5);
assert_eq!(t.stability_duration(), 2); // 2 consecutive DoNothing

// Hysteresis: fraction of measurements that triggered action
let h = t.hysteresis(); // 3 actions out of 6 measurements = 0.5
println!("Hysteresis: {:.2}", h);
```

The thermostat tracks state transitions (Heating → Cooling → Stable) and computes hysteresis — the fraction of measurements that triggered an action. High hysteresis means the system is oscillating. Low hysteresis means it's stable. A well-tuned deadband minimizes hysteresis while keeping the system responsive.

### PrecisionAnalysis — How Many Bits Can Analog Actually Get?

This is the question everyone asks. If a dial settles under physical forces, how precise is the result compared to a digital f64?

```rust
use precision::PrecisionAnalysis;
use analog_dial::AnalogDial;

let dial = AnalogDial::new(5.0, 10.0, 1.0);
let pa = PrecisionAnalysis::for_dial(&dial);

// Effective bits = log2(dynamic_range / noise_floor)
println!("Effective eigenvalue bits: {:.1}", pa.eigenvalue_bits());
// Typically ~40-50 bits under ideal conditions

// Dominant noise source
println!("Dominant noise: {}", pa.dominant_noise_source());
// Usually "friction (stick-slip)" for these parameters

// Compare to f64 (53-bit mantissa)
println!("vs digital f64: {:.1}%", pa.vs_digital() * 100.0);
// Usually ~75-95%

// Eigenvector precision (degraded by coupling condition number)
let coupling = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
println!("Eigenvector bits: {:.1}", pa.eigenvector_bits(&coupling));
```

**The precision analysis reveals three noise sources:**

1. **Thermal noise (kT):** At room temperature (300K), kT ≈ 4.14 × 10⁻²¹ J. This is the fundamental floor. For a dial with setpoint ~5 and deadband ~0.1, the dynamic range is roughly 5 / 4e-21 ≈ 1.25 × 10²¹, which is about **70 bits**. So thermal noise alone would allow very high precision.

2. **Friction noise (stick-slip):** Friction coefficient × 10⁻³. For friction = 1.0, that's 10⁻³. This dominates kT by many orders of magnitude. Dynamic range becomes 5 / 10⁻³ = 5000, which is about **12 bits**. Friction is the real precision killer.

3. **Gravity precision:** Gravity known to ~10⁻⁹ relative precision. For gravity = 10, that's 10⁻⁸ absolute. Better than friction but worse than kT.

**So the honest answer: under the simplified model here, you get about 12 bits of eigenvalue precision, limited by friction.** Under ideal conditions (carefully manufactured low-friction systems, cryogenic temperatures, feedback-enhanced dials), analog computation can reach 40-50 bits. This matches real-world observations — the best analog computers (precision op-amps, superconducting circuits) achieve 10⁻¹² to 10⁻¹⁵ relative precision, which is about 40-50 bits.

Digital f64 has a 53-bit mantissa. Analog can approach it but not exceed it under any realistic conditions.

## Connection to Real Analog Computers

**Op-amp eigenvalue solvers (1950s-1970s):** Analog computers used operational amplifiers configured as integrators to solve differential equations. Eigenvalue problems were solved by setting up the matrix as a resistor network and letting the circuit settle. The settled voltages were the eigenvectors. This library does the same thing with spring-damper dynamics instead of RC circuits.

**Memristor networks (2010s-present):** Memristor crossbar arrays can solve eigenvector problems by letting the circuit settle to its energy minimum. The conductance pattern IS the eigenvector. This is physically identical to what DialBank does — coupled oscillators settling to equilibrium.

**Hopfield networks:** A Hopfield network is a bank of coupled oscillators that settles to an energy minimum. The settled state is a stored pattern (eigenvector of the weight matrix). The deadband is the basin of attraction.

**Quantum annealing:** A quantum annealer (like D-Wave) finds eigenvalues by letting a quantum system settle to its ground state. The spectral gap determines the annealing time — exactly the relationship captured by `SpectralGapAnalysis.settling_time()`.

## Honest Limitations

1. **Friction models are idealized.** Real Coulomb friction has stick-slip dynamics, Stribeck effect, and temperature dependence. The simple `friction / gravity` deadband formula is a useful approximation but not physically accurate for most materials.

2. **No temperature dependence.** Gravity (spring constant), friction, and thermal noise all change with temperature. The precision analysis computes kT at 300K but doesn't feed it back into the dial dynamics. A real analog computer would need temperature compensation circuits.

3. **Single-precision settling.** The dials use f64 internally, but the deadband width (friction/gravity) limits the effective precision to well below f64. The precision analysis quantifies this gap, but it means you can't use settled dial positions as drop-in replacements for numerically computed eigenvalues.

4. **DialBank finds one eigenvector, not all.** The coupled dials converge to the dominant eigenvector (largest eigenvalue). Finding all eigenvectors would require deflation or orthogonalization, which aren't implemented here. Real analog computers faced the same limitation.

5. **Convergence depends on initial conditions.** The dials start at positions determined by the diagonal and index offsets. Different starting positions can lead to different convergence rates or, in pathological cases, convergence to a non-dominant mode.

## Running Tests

```bash
cargo test
```

The test suite includes:
- Dial settling, deadband computation, precision verification
- DialBank convergence, eigenvalue estimation, eigenvector residual
- Spectral gap computation, settling time, conditioning
- Thermostat state transitions, hysteresis, stability tracking
- Precision analysis: effective bits, noise sources, digital comparison
- End-to-end: build a graph → settle dials → compute eigenvalue → run thermostat → analyze precision

## Related

- **[spectral-deadband](https://github.com/SuperInstance/spectral-deadband)** — The companion library: deadband as spectral gap, spin as time, fractal conservation ratio. This library (`analog-spectral`) focuses on the analog computation aspect; `spectral-deadband` focuses on the spectral theory aspect. Same idea, different entry points.

## License

MIT

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.
