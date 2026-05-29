/// A physical dial under gravity.
/// Position = eigenvalue estimate
/// Deadband = spectral gap (gravity/friction equilibrium)
/// Settling time = convergence rate
pub struct AnalogDial {
    pub position: f64,
    pub velocity: f64,
    pub setpoint: f64,
    pub gravity: f64,
    pub friction: f64,
    pub mass: f64,
}

impl AnalogDial {
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

    /// One timestep of analog dynamics.
    /// Within deadband: friction > gravity component → no movement.
    /// Outside deadband: gravity pulls toward setpoint.
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

    /// Settle to equilibrium. Returns number of steps to converge.
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

    /// Deadband width = friction / gravity (physical spectral gap).
    pub fn deadband(&self) -> f64 {
        if self.gravity > 0.0 {
            self.friction / self.gravity
        } else {
            f64::INFINITY
        }
    }

    /// Is the dial settled (within deadband)?
    pub fn is_settled(&self) -> bool {
        (self.position - self.setpoint).abs() <= self.deadband() && self.velocity.abs() < self.deadband()
    }

    /// Precision of the settled value (related to deadband width).
    pub fn precision(&self) -> f64 {
        self.deadband()
    }
}
