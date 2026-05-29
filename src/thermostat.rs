/// A spectral thermostat — uses deadband logic to control a system.
/// This is what a real thermostat does, but formalized.

#[derive(Debug, Clone, PartialEq)]
pub enum ThermostatState {
    Heating,
    Cooling,
    Stable,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Action {
    IncreaseCR,
    DecreaseCR,
    DoNothing,
}

pub struct SpectralThermostat {
    pub setpoint: f64,
    pub deadband: f64,
    pub current: f64,
    pub history: Vec<f64>,
    pub action_history: Vec<Action>,
    state: ThermostatState,
    stable_count: usize,
}

impl SpectralThermostat {
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

    /// Measure and decide action based on deadband logic.
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

    /// Steps since last non-DoNothing action.
    pub fn stability_duration(&self) -> usize {
        self.stable_count
    }

    /// Hysteresis: fraction of measurements that triggered an action.
    pub fn hysteresis(&self) -> f64 {
        if self.action_history.is_empty() {
            return 0.0;
        }
        let actions = self.action_history.iter()
            .filter(|a| !matches!(a, Action::DoNothing))
            .count();
        actions as f64 / self.action_history.len() as f64
    }

    /// Current state.
    pub fn state(&self) -> &ThermostatState {
        &self.state
    }
}
