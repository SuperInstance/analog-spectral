//! # analog-spectral
//!
//! **Eigenvalue estimation as mechanical settling — physical dials converge under
//! gravity, friction creates deadbands equal to spectral gaps.**
//!
//! This library treats eigenvalue computation as a physical process. Each dial
//! represents an eigenvalue estimate that converges to its setpoint (true eigenvalue)
//! under spring-like restoring forces (gravity), with friction creating deadbands
//! that equal the spectral gap between eigenvalues.
//!
//! # Key Insight
//!
//! A damped harmonic oscillator settling to equilibrium IS eigenvalue iteration.
//! The deadband (friction/gravity) is the spectral gap — the region where the
//! eigenvalue estimate is "close enough" that further convergence is noise-dominated.
//!
//! # Modules
//!
//! - [`analog_dial`] — Single dial under gravity: the atomic unit of spectral computation
//! - [`dial_bank`] — N coupled dials for eigenvector/eigenvalue estimation
//! - [`spectral_gap`] — Eigenvalue gaps as deadband widths
//! - [`thermostat`] — Deadband-based control for spectral computation
//! - [`precision`] — Analog precision analysis vs digital limits

pub mod analog_dial;
pub mod dial_bank;
pub mod spectral_gap;
pub mod thermostat;
pub mod precision;

pub use analog_dial::AnalogDial;
pub use dial_bank::DialBank;
pub use spectral_gap::SpectralGapAnalysis;
pub use thermostat::{SpectralThermostat, Action, ThermostatState};
pub use precision::PrecisionAnalysis;
