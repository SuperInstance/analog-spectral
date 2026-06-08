# analog-spectral

> Analog eigenvalue computation — dials settle under gravity to spectral decomposition

[![crates.io](https://img.shields.io/crates/v/analog-spectral.svg)](https://crates.io/crates/analog-spectral)
[![docs.rs](https://docs.rs/analog-spectral/badge.svg)](https://docs.rs/analog-spectral)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## The Problem

Eigenvalue decomposition is fundamental to linear algebra — it reveals the principal directions and magnitudes of a linear transformation. Standard algorithms (QR iteration, Lanczos, power iteration) treat this as a purely algebraic problem. But eigenvalues have a natural physical interpretation: they're the equilibrium frequencies of a coupled oscillator system.

What if we modeled spectral decomposition as a physical simulation?

## The Insight

A symmetric matrix **A** can be thought of as describing the coupling between oscillating weights ("dials") connected by springs. When you place these dials at arbitrary positions and let physics take over:

1. The coupling forces (from **A** × **x**) push dials toward alignment with eigenvectors
2. Damping dissipates kinetic energy, preventing oscillation
3. Gravity (the Rayleigh quotient) provides the eigenvalue estimate
4. The system settles — eigenvalues emerge from the equilibrium positions

The **spectral gap** (distance between eigenvalues) directly determines how fast the dials settle. A large gap means fast convergence. A degenerate matrix (repeated eigenvalues) means slow, ambiguous settling — just like in the real physics.

## How It Works

### The Physical Model

Each **Dial** is a mass on a spring with position, velocity, mass, and damping:

```
Force = A·x - λ·x    (residual force from matrix coupling)
Damping = -c·v         (energy dissipation)
Acceleration = Force / mass
```

### The Algorithm

1. **Initialize** dials at arbitrary positions (like random starting angles on a physical dial)
2. **Iterate** each time step:
   - Compute Rayleigh quotient λ = (xᵀAx)/(xᵀx) — the "gravity" reading
   - Compute forces F = Ax - λx — residual drives dials toward eigenvectors
   - Apply forces with damping — physics integration step
   - Normalize positions — prevents numerical overflow (like power iteration)
3. **Converge** when residual ‖Ax - λx‖ falls below tolerance
4. **Deflate** — subtract found eigenpair, repeat for next eigenvalue

### Why It Works

This is power iteration in disguise — but dressed as physics. The normalization step is exactly the renormalization in power iteration. The damping ensures convergence rather than oscillation. The Rayleigh quotient gives the optimal eigenvalue estimate at each step.

The physical framing isn't just aesthetic: it suggests natural extensions (variable damping schedules, non-linear springs for non-symmetric matrices, thermal noise for escaping local minima).

## Code

```rust
use analog_spectral::prelude::*;

// Define a symmetric matrix
let field = GravityField::new(vec![
    vec![2.0, 1.0],
    vec![1.0, 3.0],
]);

// Find the dominant eigenvalue
let settler = Settler::new(SettleConfig::default());
let system = DialSystem::new_random(2, 1.0);
let result = settler.settle_dominant(&field, system);

println!("Dominant eigenvalue: {:?}", result.dominant_eigenvalue());
// → ~3.618 (the golden ratio + 2)

// Find all eigenvalues via deflation
let all = settler.settle_all(&field, 2);
println!("All eigenvalues: {:?}", all.eigenvalues);
// → [~3.618, ~1.382]
```

Or use the high-level solver:

```rust
let solver = SpectralSolver::new(matrix)
    .with_tolerance(1e-8)
    .with_max_iterations(5000);

let result = solver.solve(3); // Find top 3 eigenvalues
```

## Module Map

| Module | Purpose |
|---|---|
| `dial` | `Dial` struct — position, velocity, mass, damping. Physics integration (semi-implicit Euler) |
| `system` | `DialSystem` — collection of coupled dials with global damping and coupling strength |
| `gravity` | `GravityField` — the matrix encoded as coupling forces. Matrix-vector products, Rayleigh quotient, residual computation |
| `settle` | `Settler` + `SettleConfig` — the settling algorithm. Iterates dial dynamics until convergence |
| `spectral` | `SpectralSolver` + `SpectralResult` — high-level API. Dominant eigenvalue, multi-eigenvalue via deflation, result verification |
| `convergence` | `ConvergenceTracker` + `ConvergenceInfo` — tracks residual, convergence rate, spectral gap estimation |

## Design Decisions

**Pure Rust, no external math libs.** Matrix operations (matvec, Rayleigh quotient, residual) are implemented inline. For a small spectral library, pulling in `nalgebra` or `ndarray` is overkill — the core operations are ~20 lines each.

**serde + thiserror.** The only dependencies. `serde` for serializing dial states and results (useful for checkpointing long simulations). `thiserror` for clean error types.

**Power iteration semantics.** The algorithm naturally implements power iteration with Rayleigh quotient refinement. This means it finds the dominant eigenvalue first and best. Smaller eigenvalues found via deflation will be less precise — this is a known trade-off.

**Semi-implicit Euler integration.** Simple, stable, and energy-preserving-enough for our purposes. The damping ensures we don't need a more sophisticated integrator.

**No `unsafe`.** The entire crate is safe Rust. Performance-critical code could benefit from SIMD or BLAS calls, but for the educational/prototyping use case, clarity wins.

**Configurable physics.** Mass, damping, coupling strength, and time step are all tunable. This lets you experiment with convergence behavior — heavier dials settle slower, more damping prevents overshooting.

## Limitations

- **Symmetric matrices work best.** Non-symmetric matrices can have complex eigenvalues, which this real-valued simulation can't capture directly.
- **Deflation accumulates error.** Each successive eigenvalue is slightly less accurate than the previous one.
- **Negative dominant eigenvalues.** Power iteration oscillates when the dominant eigenvalue is negative. The Rayleigh quotient approach may converge to a different eigenvalue in this case.
- **Not a replacement for LAPACK.** This is an educational/experimental library. For production eigenvalue computations, use `nalgebra` or `rust-ndarray` with LAPACK bindings.

## License

MIT
