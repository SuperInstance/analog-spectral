//! Analog eigenvalue computation via physical dials

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
