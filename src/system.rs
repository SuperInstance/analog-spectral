//! DialSystem: a collection of coupled dials forming an analog spectral system.

use crate::dial::Dial;
use serde::{Deserialize, Serialize};

/// A system of coupled dials representing the spectral decomposition problem.
///
/// The dial system maintains N dials whose positions converge to reveal
/// eigenvectors and eigenvalues of the coupling matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DialSystem {
    /// The dials in this system.
    dials: Vec<Dial>,
    /// Global damping applied to all dials.
    global_damping: f64,
    /// Coupling strength between dials.
    coupling_strength: f64,
}

impl DialSystem {
    /// Create a new dial system with n dials initialized to given positions.
    pub fn new(positions: Vec<f64>) -> Self {
        let dials = positions.into_iter().map(Dial::new).collect();
        Self {
            dials,
            global_damping: 0.5,
            coupling_strength: 1.0,
        }
    }

    /// Create a system with n dials at random-ish initial positions.
    pub fn new_random(n: usize, seed: f64) -> Self {
        let positions: Vec<f64> = (0..n)
            .map(|i| {
                // Simple deterministic spread, not cryptographic randomness
                let t = seed + i as f64;
                (t * 2.39996).sin() * 5.0 + (t * 1.71828).cos() * 3.0
            })
            .collect();
        Self::new(positions)
    }

    /// Create a system with n dials at origin.
    pub fn new_zeros(n: usize) -> Self {
        Self::new(vec![0.0; n])
    }

    /// Set global damping.
    pub fn with_global_damping(mut self, damping: f64) -> Self {
        self.global_damping = damping;
        self
    }

    /// Set coupling strength.
    pub fn with_coupling_strength(mut self, strength: f64) -> Self {
        self.coupling_strength = strength;
        self
    }

    /// Get a reference to the dials.
    pub fn dials(&self) -> &[Dial] {
        &self.dials
    }

    /// Get a mutable reference to the dials.
    pub fn dials_mut(&mut self) -> &mut Vec<Dial> {
        &mut self.dials
    }

    /// Number of dials.
    pub fn len(&self) -> usize {
        self.dials.len()
    }

    /// Whether the system has no dials.
    pub fn is_empty(&self) -> bool {
        self.dials.is_empty()
    }

    /// Get the current positions of all dials as a vector.
    pub fn positions(&self) -> Vec<f64> {
        self.dials.iter().map(|d| d.position).collect()
    }

    /// Get the current velocities of all dials.
    pub fn velocities(&self) -> Vec<f64> {
        self.dials.iter().map(|d| d.velocity).collect()
    }

    /// Total kinetic energy of the system.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.dials.iter().map(|d| d.kinetic_energy()).sum()
    }

    /// Check if all dials have settled (velocity below threshold).
    pub fn all_settled(&self, velocity_threshold: f64) -> bool {
        self.dials.iter().all(|d| d.is_settled(velocity_threshold))
    }

    /// Normalize dial positions to unit vector.
    pub fn normalize_positions(&mut self) {
        let norm: f64 = self
            .dials
            .iter()
            .map(|d| d.position * d.position)
            .sum::<f64>()
            .sqrt();
        if norm > 1e-30 {
            for dial in &mut self.dials {
                dial.position /= norm;
            }
        }
    }

    /// Apply forces to all dials for one time step.
    pub fn apply_forces(&mut self, forces: &[f64], dt: f64) {
        assert_eq!(forces.len(), self.dials.len());
        for (dial, &force) in self.dials.iter_mut().zip(forces.iter()) {
            let total_force = self.coupling_strength * force;
            dial.step_coupled(total_force, dt);
        }
    }

    /// Apply damping to all dials.
    pub fn apply_damping(&mut self) {
        for dial in &mut self.dials {
            dial.velocity *= 1.0 - self.global_damping * 0.01;
        }
    }

    /// Reset all dials to given positions.
    pub fn reset(&mut self, positions: Vec<f64>) {
        assert_eq!(positions.len(), self.dials.len());
        for (dial, pos) in self.dials.iter_mut().zip(positions) {
            dial.reset(pos);
        }
    }

    /// Perturb dial positions slightly to break symmetry.
    pub fn perturb(&mut self, amount: f64) {
        for (i, dial) in self.dials.iter_mut().enumerate() {
            let offset = amount * ((i as f64 + 1.0) * 0.1).sin();
            dial.position += offset;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_system_correct_count() {
        let sys = DialSystem::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(sys.len(), 3);
        assert!(!sys.is_empty());
    }

    #[test]
    fn positions_match_init() {
        let sys = DialSystem::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(sys.positions(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn velocities_initially_zero() {
        let sys = DialSystem::new(vec![1.0, 2.0]);
        assert_eq!(sys.velocities(), vec![0.0, 0.0]);
    }

    #[test]
    fn total_kinetic_initially_zero() {
        let sys = DialSystem::new(vec![5.0, 5.0]);
        assert_eq!(sys.total_kinetic_energy(), 0.0);
    }

    #[test]
    fn all_settled_initially() {
        let sys = DialSystem::new(vec![1.0, 2.0]);
        assert!(sys.all_settled(0.1));
    }

    #[test]
    fn normalize_positions() {
        let mut sys = DialSystem::new(vec![3.0, 4.0]);
        sys.normalize_positions();
        let pos = sys.positions();
        let norm = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn apply_forces_moves_dials() {
        let mut sys = DialSystem::new(vec![0.0, 0.0]);
        sys.apply_forces(&[1.0, 1.0], 0.1);
        // Should have moved
        assert!(sys.positions()[0] > 0.0);
    }

    #[test]
    fn reset_clears_state() {
        let mut sys = DialSystem::new(vec![1.0, 2.0]);
        sys.apply_forces(&[10.0, 10.0], 0.1);
        sys.reset(vec![0.0, 0.0]);
        assert_eq!(sys.positions(), vec![0.0, 0.0]);
        assert_eq!(sys.velocities(), vec![0.0, 0.0]);
    }

    #[test]
    fn perturb_adds_offset() {
        let mut sys = DialSystem::new(vec![1.0, 1.0]);
        let before = sys.positions();
        sys.perturb(0.1);
        let after = sys.positions();
        // Positions should have changed
        assert_ne!(before, after);
    }

    #[test]
    fn zeros_system() {
        let sys = DialSystem::new_zeros(4);
        assert_eq!(sys.len(), 4);
        assert_eq!(sys.positions(), vec![0.0; 4]);
    }

    #[test]
    fn empty_system() {
        let sys = DialSystem::new(vec![]);
        assert!(sys.is_empty());
    }
}
