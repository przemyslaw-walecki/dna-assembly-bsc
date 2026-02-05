use std::fs;
use std::path::{Path, PathBuf};

pub fn make_tmp_dir(tag: &str) -> PathBuf {
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();

    let mut dir = std::env::temp_dir();
    dir.push(format!("rust_bubble_labelling_it_{}_{}_{}", tag, pid, nanos));
    fs::create_dir_all(&dir).unwrap();
    dir
}

pub fn write_file(dir: &Path, name: &str, content: &str) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, content).unwrap();
    p
}

pub fn read_lines(path: &Path) -> Vec<String> {
    let s = fs::read_to_string(path).unwrap_or_default();
    s.lines().map(|l| l.to_string()).collect()
}

/// Uruchamia binarkę z argumentami; zwraca (status_code, stdout, stderr).
pub fn run_bin(args: &[&str]) -> (i32, String, String) {
    // Jeśli binarka nazywa się inaczej niż pakiet, zmień suffix env var.
    let exe = env!("CARGO_BIN_EXE_rust_bubble_labelling");

    let out = std::process::Command::new(exe)
        .args(args)
        // Ułatwia testy w CI/terminalach bez TTY (indicatif mniej „szaleje”).
        .env("TERM", "dumb")
        .output()
        .expect("failed to run binary");

    let code = out.status.code().unwrap_or(-1);
    let stdout = String::from_utf8_lossy(&out.stdout).to_string();
    let stderr = String::from_utf8_lossy(&out.stderr).to_string();

    (code, stdout, stderr)
}
