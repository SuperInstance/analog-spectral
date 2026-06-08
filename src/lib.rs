//! # analog-spectral
//!
//! Analog eigenvalue computation where dials settle under gravity.
//!
//! Models spectral analysis as a physical system where "dials" (oscillating
//! weights on springs) settle to eigenvalues under gravity. The spectral gap
//! determines convergence rate.

mod convergence;
mod dial;
mod gravity;
mod settle;
mod spectral;
mod system;

pub use convergence::{ConvergenceInfo, ConvergenceTracker};
pub use dial::Dial;
pub use gravity::GravityField;
pub use settle::{SettleConfig, Settler};
pub use spectral::{SpectralResult, SpectralSolver};
pub use system::DialSystem;

/// Re-export core types for convenience.
pub mod prelude {
    pub use crate::{
        ConvergenceInfo, ConvergenceTracker, Dial, DialSystem, GravityField, SettleConfig,
        Settler, SpectralResult, SpectralSolver,
    };
}
