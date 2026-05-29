fn main() {
    println!("analog-spectral: physical dials, gravity deadbands, spectral computation");
}

mod analog_dial;
mod dial_bank;
mod spectral_gap;
mod thermostat;
mod precision;

#[cfg(test)]
mod tests {
    use super::*;
    use analog_dial::AnalogDial;
    use dial_bank::DialBank;
    use spectral_gap::SpectralGapAnalysis;
    use thermostat::{SpectralThermostat, Action, ThermostatState};
    use precision::PrecisionAnalysis;

    // === AnalogDial tests ===

    #[test]
    fn dial_settles_to_setpoint_within_deadband() {
        let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
        dial.position = 0.0;
        dial.settle(0.01, 1e-10);
        assert!(dial.is_settled());
        assert!((dial.position - 5.0).abs() < dial.deadband());
    }

    #[test]
    fn dial_deadband_equals_friction_over_gravity() {
        let dial = AnalogDial::new(3.0, 8.0, 2.0);
        let expected = 2.0 / 8.0;
        assert!((dial.deadband() - expected).abs() < 1e-12);
    }

    #[test]
    fn dial_settling_time_inversely_proportional_to_gravity() {
        // Measure settling in simulated time (steps × dt), not just steps
        let dt = 0.001;
        let tolerance = 0.1;

        let mut d1 = AnalogDial::new(10.0, 4.0, 0.02);
        d1.position = 0.0;
        d1.velocity = 0.0;
        let steps1 = d1.settle(dt, tolerance);
        let time1 = steps1 as f64 * dt;

        let mut d2 = AnalogDial::new(10.0, 16.0, 0.02);
        d2.position = 0.0;
        d2.velocity = 0.0;
        let steps2 = d2.settle(dt, tolerance);
        let time2 = steps2 as f64 * dt;

        // Higher gravity → faster settling (less simulated time)
        assert!(time2 < time1, "d1 time={}, d2 time={}", time1, time2);
    }

    #[test]
    fn dial_precision_related_to_deadband() {
        let dial = AnalogDial::new(5.0, 10.0, 1.0);
        let precision = dial.precision();
        assert!(precision > 0.0);
        assert!((precision - dial.deadband()).abs() < 1e-12);
    }

    #[test]
    fn dial_step_moves_toward_setpoint() {
        let mut dial = AnalogDial::new(10.0, 5.0, 0.1);
        dial.position = 0.0;
        dial.velocity = 0.0;
        let prev = dial.position;
        dial.step(0.1);
        // Should move toward setpoint (positive direction)
        assert!(dial.position > prev || dial.velocity > 0.0);
    }

    // === DialBank tests ===

    #[test]
    fn dial_bank_coupled_dials_converge() {
        // Simple 2x2 symmetric matrix: [[2, 1], [1, 2]]
        let coupling = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let mut bank = DialBank::new(2, coupling.clone());
        bank.settle(0.01, 1e-8);
        let positions = bank.read_positions();
        // Should converge to something (not diverge)
        for p in &positions {
            assert!(p.is_finite());
        }
    }

    #[test]
    fn dial_bank_eigenvalue_estimate() {
        // 2x2 identity: eigenvalues are both 1.0
        let coupling = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let mut bank = DialBank::new(2, coupling.clone());
        bank.settle(0.01, 1e-8);
        let ev = bank.eigenvalue_estimate();
        assert!(ev.is_finite());
    }

    #[test]
    fn dial_bank_verify_eigenvector() {
        // Symmetric matrix [[3, 1], [1, 3]] — eigenvalues 4, 2
        let adj = vec![vec![3.0, 1.0], vec![1.0, 3.0]];
        let mut bank = DialBank::new(2, adj.clone());
        bank.settle(0.01, 1e-8);
        let residual = bank.verify_eigenvector(&adj);
        assert!(residual < 0.5, "Residual too large: {}", residual);
    }

    #[test]
    fn dial_bank_step_simultaneous() {
        let coupling = vec![vec![2.0, 0.0], vec![0.0, 2.0]];
        let mut bank = DialBank::new(2, coupling);
        bank.step(0.01);
        // Just check it doesn't crash and positions are finite
        for p in bank.read_positions() {
            assert!(p.is_finite());
        }
    }

    // === SpectralGap tests ===

    #[test]
    fn spectral_gaps_computed_correctly() {
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);
        assert!((analysis.gaps[0] - 1.0).abs() < 1e-12);
        assert!((analysis.gaps[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn spectral_largest_gap() {
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);
        let (idx, gap) = analysis.largest_gap();
        assert_eq!(idx, 1);
        assert!((gap - 3.0).abs() < 1e-12);
    }

    #[test]
    fn spectral_settling_time_inversely_proportional() {
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);
        // Larger gap → shorter settling time
        let t0 = analysis.settling_time(0);
        let t1 = analysis.settling_time(1);
        assert!(t1 < t0, "Gap 1 (3.0) should settle faster than gap 0 (1.0)");
    }

    #[test]
    fn spectral_conditioning() {
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 2.0, 5.0]);
        let cond = analysis.conditioning();
        // min_gap/max_gap = 1.0/3.0
        assert!((cond - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn spectral_deadband_at() {
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![1.0, 3.0, 7.0]);
        let db = analysis.deadband_at(0);
        assert!(db > 0.0);
        // Should be gap / max_eigenvalue
        let expected = 2.0 / 7.0;
        assert!((db - expected).abs() < 1e-12);
    }

    // === Thermostat tests ===

    #[test]
    fn thermostat_do_nothing_within_deadband() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        let action = t.measure(0.5);
        assert!(matches!(action, Action::DoNothing));
    }

    #[test]
    fn thermostat_increase_cr_below_setpoint() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        let action = t.measure(0.3);
        assert!(matches!(action, Action::IncreaseCR));
    }

    #[test]
    fn thermostat_decrease_cr_above_setpoint() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        let action = t.measure(0.7);
        assert!(matches!(action, Action::DecreaseCR));
    }

    #[test]
    fn thermostat_stability_duration_increases() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        t.measure(0.5);
        assert_eq!(t.stability_duration(), 1);
        t.measure(0.5);
        assert_eq!(t.stability_duration(), 2);
    }

    #[test]
    fn thermostat_hysteresis() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        t.measure(0.5);
        t.measure(0.7);
        let h = t.hysteresis();
        assert!(h >= 0.0 && h <= 1.0);
    }

    #[test]
    fn thermostat_state_transitions() {
        let mut t = SpectralThermostat::new(0.5, 0.1);
        t.measure(0.3);
        assert!(matches!(t.state(), ThermostatState::Heating));
        t.measure(0.7);
        assert!(matches!(t.state(), ThermostatState::Cooling));
        t.measure(0.5);
        assert!(matches!(t.state(), ThermostatState::Stable));
    }

    // === PrecisionAnalysis tests ===

    #[test]
    fn precision_effective_bits() {
        let dial = AnalogDial::new(5.0, 10.0, 1.0);
        let pa = PrecisionAnalysis::for_dial(&dial);
        let bits = pa.eigenvalue_bits();
        assert!(bits > 0.0);
    }

    #[test]
    fn precision_vs_digital_less_than_one() {
        let dial = AnalogDial::new(5.0, 10.0, 1.0);
        let pa = PrecisionAnalysis::for_dial(&dial);
        let ratio = pa.vs_digital();
        assert!(ratio < 1.0, "Analog should be less precise than f64");
    }

    #[test]
    fn precision_dominant_noise_source() {
        let dial = AnalogDial::new(5.0, 10.0, 1.0);
        let pa = PrecisionAnalysis::for_dial(&dial);
        let source = pa.dominant_noise_source();
        assert!(!source.is_empty());
    }

    #[test]
    fn precision_eigenvector_bits() {
        let dial = AnalogDial::new(5.0, 10.0, 1.0);
        let pa = PrecisionAnalysis::for_dial(&dial);
        let coupling = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let bits = pa.eigenvector_bits(&coupling);
        assert!(bits > 0.0);
    }

    // === End-to-end test ===

    #[test]
    fn end_to_build_graph_compute_eigenvalues() {
        // Build a simple graph adjacency matrix (triangle + isolated node)
        // [[0, 1, 1, 0], [1, 0, 1, 0], [1, 1, 0, 0], [0, 0, 0, 0]]
        let adj = vec![
            vec![0.0, 1.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0, 0.0],
            vec![1.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
        ];

        let mut bank = DialBank::new(4, adj.clone());
        let steps = bank.settle(0.01, 1e-8);
        let positions = bank.read_positions();

        // All positions should be finite
        for p in &positions {
            assert!(p.is_finite());
        }

        // Eigenvalue estimate should be finite
        let ev = bank.eigenvalue_estimate();
        assert!(ev.is_finite());

        // Spectral gap analysis from the eigenvalue estimate
        let analysis = SpectralGapAnalysis::from_eigenvalues(vec![ev]);
        let (idx, gap) = analysis.largest_gap();
        assert_eq!(idx, 0);
        assert!(gap >= 0.0);

        // Run thermostat with the result
        let mut thermostat = SpectralThermostat::new(2.0, 0.5);
        let action = thermostat.measure(ev);
        // Whatever the action, it should be valid
        match action {
            Action::IncreaseCR | Action::DecreaseCR | Action::DoNothing => {}
        }

        // Precision analysis
        let pa = PrecisionAnalysis::for_dial(&AnalogDial::new(ev, 10.0, 1.0));
        assert!(pa.vs_digital() < 1.0);
    }
}
