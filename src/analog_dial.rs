//! Single analog dial under gravity — the atomic unit of spectral computation.
//!
//! A dial's position represents an eigenvalue estimate. Gravity provides a
//! restoring force toward the setpoint (true eigenvalue), while friction
//! creates a deadband equal to the spectral gap. The settling dynamics
//! mirror iterative eigenvalue convergence.

/// A physical dial under gravity.
///
/// The dial's position converges to its setpoint under spring-like restoring
/// forces. The deadband (friction/gravity) acts as a natural spectral gap —
/// within it, the dial is "close enough" and friction dominates.
///
/// # Physical Analogy
///
/// - **Position** = eigenvalue estimate
/// - **Deadband** = spectral gap (gravity/friction equilibrium)
/// - **Settling time** = convergence rate of the eigenvalue algorithm
pub struct AnalogDial {
    /// Current position — the eigenvalue estimate.
    pub position: f64,
    /// Current velocity — rate of change of the estimate.
    pub velocity: f64,
    /// Target position — the true eigenvalue being sought.
    pub setpoint: f64,
    /// Restoring force strength — higher gravity = faster convergence.
    pub gravity: f64,
    /// Damping coefficient — creates deadband (spectral gap = friction/gravity).
    pub friction: f64,
    /// Effective mass — controls inertia of the dial.
    pub mass: f64,
}

impl AnalogDial {
    /// Create a new dial initialized at its setpoint with zero velocity.
    ///
    /// # Arguments
    ///
    /// * `setpoint` — Target eigenvalue the dial seeks
    /// * `gravity` — Restoring force strength (higher → faster convergence)
    /// * `friction` — Damping coefficient (deadband = friction/gravity)
    pub fn new(setpoint: f64, gravity: f64, friction: f64) -> AnalogDial {
        AnalogDial {
            position: setpoint,
            velocity: 0.0,
            setpoint,
            gravity,
            friction,
            mass: 1.0,
        }
    }

    /// Advance one timestep of analog dynamics.
    ///
    /// Within the deadband: friction dominates and the dial stops moving.
    /// Outside the deadband: gravity pulls the dial toward the setpoint,
    /// damped by friction opposing velocity.
    pub fn step(&mut self, dt: f64) {
        let displacement = self.position - self.setpoint;

        // Restoring force (gravity-like, proportional to displacement)
        let spring_force = -self.gravity * displacement;

        // Damping (friction) opposes velocity
        let damping_force = -self.friction * self.velocity;

        // Coulomb-like friction: if nearly stopped and within deadband, stop
        let deadband = self.deadband();
        if displacement.abs() < deadband && self.velocity.abs() < deadband {
            self.velocity = 0.0;
            return;
        }

        let acceleration = (spring_force + damping_force) / self.mass;
        self.velocity += acceleration * dt;
        self.position += self.velocity * dt;
    }

    /// Settle the dial to equilibrium by iterating until convergence.
    ///
    /// Returns the number of steps needed to converge within `tolerance`.
    /// Aborts after 10 million steps if convergence is not achieved.
    pub fn settle(&mut self, dt: f64, tolerance: f64) -> usize {
        let max_steps = 10_000_000;
        for i in 0..max_steps {
            let prev = self.position;
            self.step(dt);

            // Check convergence
            let change = (self.position - prev).abs();
            let offset = (self.position - self.setpoint).abs();

            if change < tolerance && offset < tolerance {
                return i + 1;
            }
        }
        max_steps
    }

    /// Compute the deadband width = friction / gravity.
    ///
    /// This is the physical analog of the spectral gap — the region where
    /// the eigenvalue estimate is "close enough" that numerical noise
    /// dominates any further convergence.
    pub fn deadband(&self) -> f64 {
        if self.gravity > 0.0 {
            self.friction / self.gravity
        } else {
            f64::INFINITY
        }
    }

    /// Check if the dial is settled within its deadband.
    pub fn is_settled(&self) -> bool {
        (self.position - self.setpoint).abs() <= self.deadband() && self.velocity.abs() < self.deadband()
    }

    /// Precision of the settled value (equal to deadband width).
    ///
    /// Smaller deadband = higher precision eigenvalue estimate,
    /// but requires lower friction or higher gravity.
    pub fn precision(&self) -> f64 {
        self.deadband()
    }
}
