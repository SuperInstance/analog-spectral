//! Integration tests for analog-spectral

use analog_spectral::*;

#[test]
fn test_dial_settles() {
    let mut dial = AnalogDial::new(5.0, 10.0, 1.0);
    dial.position = 0.0;
    dial.settle(0.01, 1e-10);
    assert!(dial.is_settled());
    assert!((dial.position - 5.0).abs() < dial.deadband() + 0.1);
}

#[test]
fn test_spectral_gap_analysis() {
    let analysis = SpectralGapAnalysis::from_eigenvalues(vec![0.0, 1.0, 3.0, 6.0]);
    assert_eq!(analysis.gaps.len(), 3);
    assert_eq!(analysis.eigenvalues.len(), 4);
    let (idx, gap) = analysis.largest_gap();
    assert_eq!(idx, 2);
    assert!((gap - 3.0).abs() < 1e-10);
}

#[test]
fn test_spectral_gap_deadbands() {
    let analysis = SpectralGapAnalysis::from_eigenvalues(vec![0.0, 0.1, 5.0]);
    let db = analysis.deadband_at(0);
    assert!(db > 0.0);
    assert!(db <= 1.0);
}

#[test]
fn test_dial_bank_creation() {
    let bank = DialBank::new(3, vec![vec![0.0; 3]; 3]);
    assert_eq!(bank.dials.len(), 3);
}
