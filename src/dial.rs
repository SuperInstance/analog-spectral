//! A single dial: oscillating weight on a spring that settles under gravity.

use serde::{Deserialize, Serialize};

/// A dial representing one degree of freedom in the analog spectral system.
///
/// Each dial has a position (current value), velocity, mass, and damping.
/// Under gravitational forces from coupling with other dials, it settles
/// toward an eigenvalue.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Dial {
    /// Current position of the dial (converges toward an eigenvalue).
    pub position: f64,
    /// Current velocity of the dial.
    pub velocity: f64,
    /// Mass of the dial (affects inertia and settling rate).
    pub mass: f64,
    /// Damping coefficient (energy dissipation per unit velocity).
    pub damping: f64,
    /// Spring constant (restoring force toward equilibrium).
    pub spring: f64,
}

impl Dial {
    /// Create a new dial at the given position with default parameters.
    pub fn new(position: f64) -> Self {
        Self {
            position,
            velocity: 0.0,
            mass: 1.0,
            damping: 0.5,
            spring: 1.0,
        }
    }

    /// Create a dial with custom physical parameters.
    pub fn with_params(position: f64, mass: f64, damping: f64, spring: f64) -> Self {
        Self {
            position,
            velocity: 0.0,
            mass,
            damping,
            spring,
        }
    }

    /// Apply a force for one time step using semi-implicit Euler integration.
    ///
    /// The equation of motion is:
    /// `m * a = F_applied - damping * v - spring * (x - x_rest)`
    ///
    /// Returns the new position after the step.
    pub fn step(&mut self, force: f64, dt: f64) {
        let acceleration = (force - self.damping * self.velocity - self.spring * self.position)
            / self.mass;
        self.velocity += acceleration * dt;
        self.position += self.velocity * dt;
    }

    /// Apply force without spring restoring term (for coupled systems where
    /// the coupling matrix provides the restoring force).
    pub fn step_coupled(&mut self, force: f64, dt: f64) {
        let acceleration = (force - self.damping * self.velocity) / self.mass;
        self.velocity += acceleration * dt;
        self.position += self.velocity * dt;
    }

    /// Compute the kinetic energy of the dial.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity * self.velocity
    }

    /// Check if the dial has essentially stopped moving.
    pub fn is_settled(&self, velocity_threshold: f64) -> bool {
        self.velocity.abs() < velocity_threshold
    }

    /// Reset the dial to a new position with zero velocity.
    pub fn reset(&mut self, position: f64) {
        self.position = position;
        self.velocity = 0.0;
    }

    /// Return the current speed (absolute velocity).
    pub fn speed(&self) -> f64 {
        self.velocity.abs()
    }
}

impl Default for Dial {
    fn default() -> Self {
        Self::new(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_dial_has_zero_velocity() {
        let d = Dial::new(3.0);
        assert_eq!(d.position, 3.0);
        assert_eq!(d.velocity, 0.0);
        assert_eq!(d.mass, 1.0);
    }

    #[test]
    fn default_dial_is_at_origin() {
        let d = Dial::default();
        assert_eq!(d.position, 0.0);
    }

    #[test]
    fn step_with_zero_force_spring_settles() {
        let mut d = Dial::new(5.0);
        d.damping = 2.0;
        d.spring = 1.0;
        for _ in 0..1000 {
            d.step(0.0, 0.01);
        }
        assert!(d.position.abs() < 0.1);
        assert!(d.speed() < 0.01);
    }

    #[test]
    fn step_coupled_no_spring_restore() {
        let mut d = Dial::new(1.0);
        d.damping = 1.0;
        // With no spring, position drifts under force but damping slows it
        d.step_coupled(1.0, 0.01);
        assert!(d.position > 1.0);
    }

    #[test]
    fn kinetic_energy_zero_when_still() {
        let d = Dial::new(5.0);
        assert_eq!(d.kinetic_energy(), 0.0);
    }

    #[test]
    fn kinetic_energy_positive_when_moving() {
        let mut d = Dial::new(0.0);
        d.velocity = 3.0;
        let ke = d.kinetic_energy();
        assert!((ke - 4.5).abs() < 1e-10);
    }

    #[test]
    fn is_settled_checks_velocity() {
        let mut d = Dial::new(0.0);
        assert!(d.is_settled(0.1));
        d.velocity = 0.5;
        assert!(!d.is_settled(0.1));
    }

    #[test]
    fn reset_clears_velocity() {
        let mut d = Dial::new(5.0);
        d.velocity = 10.0;
        d.reset(2.0);
        assert_eq!(d.position, 2.0);
        assert_eq!(d.velocity, 0.0);
    }

    #[test]
    fn with_params_sets_all() {
        let d = Dial::with_params(1.0, 2.0, 3.0, 4.0);
        assert_eq!(d.mass, 2.0);
        assert_eq!(d.damping, 3.0);
        assert_eq!(d.spring, 4.0);
    }

    #[test]
    fn speed_is_absolute_velocity() {
        let mut d = Dial::new(0.0);
        d.velocity = -7.0;
        assert_eq!(d.speed(), 7.0);
    }

    #[test]
    fn heavier_mass_slower_acceleration() {
        let mut light = Dial::with_params(0.0, 1.0, 0.0, 0.0);
        let mut heavy = Dial::with_params(0.0, 10.0, 0.0, 0.0);
        light.step(10.0, 0.1);
        heavy.step(10.0, 0.1);
        assert!(light.position > heavy.position);
    }
}
