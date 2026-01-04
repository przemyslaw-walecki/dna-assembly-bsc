mod common;
use common::run_bin;

#[test]
fn cli_requires_4_args() {
    let (code, _stdout, stderr) = run_bin(&[]);
    assert_ne!(code, 0);
    // Nie spinamy się o dokładny tekst (progressbar/locale), tylko o sygnał.
    assert!(stderr.contains("Użycie:") || stderr.contains("Uzycie:"));
}
