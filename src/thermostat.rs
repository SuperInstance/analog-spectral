//! Spectral thermostat — deadband-based control for spectral computation.
//!
//! A thermostat monitors a spectral quantity (e.g., condition ratio) and
//! applies deadband logic to decide whether to increase or decrease
//! computational effort. Within the deadband, the system is "close enough"
//! and no action is taken.

/// State of the thermostat relative to its setpoint.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermostatState {
    /// Value is below the deadband — needs increase.
    Heating,
    /// Value is above the deadband — needs decrease.
    Cooling,
    /// Value is within the deadband — no action needed.
    Stable,
}

/// Control action decided by the thermostat.
#[derive(Debug, Clone, PartialEq)]
pub enum Action {
    /// Increase the condition ratio (value too low).
    IncreaseCR,
    /// Decrease the condition ratio (value too high).
    DecreaseCR,
    /// Within deadband — no action needed.
    DoNothing,
}

/// A spectral thermostat using deadband logic to control computation.
///
/// The thermostat monitors a single scalar quantity and decides whether
/// it needs adjustment based on a setpoint and deadband width. Within
/// the deadband, the system is considered stable.
///
/// # Physical Analogy
///
/// Like a real thermostat, this controller only acts when the value
/// drifts outside the comfort zone (deadband). This avoids oscillation
/// from noisy measurements near the setpoint.
pub struct SpectralThermostat {
    /// Target value for the controlled quantity.
    pub setpoint: f64,
    /// Half-width of the deadband around the setpoint.
    pub deadband: f64,
    /// Most recently measured value.
    pub current: f64,
    /// History of all measurements.
    pub history: Vec<f64>,
    /// History of all actions taken.
    pub action_history: Vec<Action>,
    state: ThermostatState,
    stable_count: usize,
}

impl SpectralThermostat {
    /// Create a new thermostat targeting `setpoint` with the given deadband.
    pub fn new(setpoint: f64, deadband: f64) -> SpectralThermostat {
        SpectralThermostat {
            setpoint,
            deadband,
            current: setpoint,
            history: Vec::new(),
            action_history: Vec::new(),
            state: ThermostatState::Stable,
            stable_count: 0,
        }
    }

    /// Measure a value and decide the control action.
    ///
    /// If the value is within [setpoint − deadband, setpoint + deadband],
    /// returns `DoNothing`. Otherwise returns `IncreaseCR` or `DecreaseCR`.
    pub fn measure(&mut self, cr: f64) -> Action {
        self.current = cr;
        self.history.push(cr);

        let deviation = cr - self.setpoint;

        let action = if deviation.abs() <= self.deadband {
            self.state = ThermostatState::Stable;
            self.stable_count += 1;
            Action::DoNothing
        } else if deviation < 0.0 {
            self.state = ThermostatState::Heating;
            self.stable_count = 0;
            Action::IncreaseCR
        } else {
            self.state = ThermostatState::Cooling;
            self.stable_count = 0;
            Action::DecreaseCR
        };

        self.action_history.push(action.clone());
        action
    }

    /// Number of consecutive stable measurements since the last action.
    pub fn stability_duration(&self) -> usize {
        self.stable_count
    }

    /// Hysteresis ratio: fraction of measurements that triggered an action.
    ///
    /// A higher hysteresis means the system oscillates more. Low hysteresis
    /// means the value stays within the deadband.
    pub fn hysteresis(&self) -> f64 {
        if self.action_history.is_empty() {
            return 0.0;
        }
        let actions = self.action_history.iter()
            .filter(|a| !matches!(a, Action::DoNothing))
            .count();
        actions as f64 / self.action_history.len() as f64
    }

    /// Current thermostat state.
    pub fn state(&self) -> &ThermostatState {
        &self.state
    }
}
