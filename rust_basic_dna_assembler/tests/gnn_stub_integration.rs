use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rust_basic_dna_assembler::KmerGraph;

fn tmp_dir(name: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let mut p = std::env::temp_dir();
    p.push(format!("ai_integration_test_{}_{}_{}", name, std::process::id(), ts));
    fs::create_dir_all(&p).unwrap();
    p
}

fn write_text(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f.flush().unwrap();
}

fn fixture(path: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push(path);
    p
}

fn find_python() -> String {
    // prefer python3, fallback python
    if Command::new("python3").arg("--version").status().is_ok() {
        "python3".to_string()
    } else {
        "python".to_string()
    }
}

#[test]
fn gnn_stub_end_to_end_generates_decisions_and_resolver_applies_them() {
    // Prawdziwy "diamond bubble":
    // s -> a -> t  (EC wysokie)
    // s -> b -> t  (EC niskie)  => powinno zostać usunięte
    let mut g = KmerGraph::new(3);

    let s = g.encode("AAA");
    let a = g.encode("AAT");
    let b = g.encode("AAC");
    let t = g.encode("TTT");

    g.edges.entry(s).or_default().extend([a, b]);
    g.edges.entry(a).or_default().insert(t);
    g.edges.entry(b).or_default().insert(t);
    g.edges.entry(t).or_default();

    g.out_degree.insert(s, 2);
    g.out_degree.insert(a, 1);
    g.out_degree.insert(b, 1);
    g.out_degree.insert(t, 0);

    g.in_degree.insert(s, 0);
    g.in_degree.insert(a, 1);
    g.in_degree.insert(b, 1);
    g.in_degree.insert(t, 2);

    // Preferujemy ścieżkę przez 'a'
    g.edge_counts.insert((s, a), 10);
    g.edge_counts.insert((a, t), 10);
    g.edge_counts.insert((s, b), 1);
    g.edge_counts.insert((b, t), 1);

    g.node_coverage.insert(s, 1);
    g.node_coverage.insert(a, 1);
    g.node_coverage.insert(b, 1);
    g.node_coverage.insert(t, 1);

    let dir = tmp_dir("gnn_stub");
    let gfa_path = dir.join("graph.gfa");
    let bubbles_path = dir.join("bubbles.jsonl");
    let decisions_path = dir.join("decisions.jsonl");

    g.write_gfa(gfa_path.to_string_lossy().as_ref()).unwrap();

    // BubbleGun-like input: ends = [s,t], inside = [a,b]
    let bubbles = format!(
        "{{\"bubbles\":[{{\"id\":1,\"ends\":[\"{}\",\"{}\"],\"inside\":[\"{}\",\"{}\"]}}]}}\n",
        s, t, a, b
    );
    write_text(&bubbles_path, &bubbles);

    let stub_path = fixture("tests/fixtures/stubs/gnn_stub.py");
    let python = find_python();

    let status = Command::new(&python)
        .arg(stub_path)
        .arg("--in")
        .arg(&bubbles_path)
        .arg("--out")
        .arg(&decisions_path)
        .arg("--gfa")
        .arg(&gfa_path)
        .arg("--ckpt")
        .arg("ignored.ckpt")
        .status()
        .expect("failed to run python gnn_stub");

    assert!(status.success(), "gnn_stub.py failed with status: {status}");

    let removed = g
        .resolve_bubbles_from_jsonl(decisions_path.to_string_lossy().as_ref())
        .unwrap();

    assert!(removed >= 1, "expected at least one edge removed, got {removed}");

    // ścieżka keep musi zostać
    assert!(g.edges.get(&s).unwrap().contains(&a));
    assert!(g.edges.get(&a).unwrap().contains(&t));
    
    // konkurencyjna gałąź musi być odcięta od startu
    assert!(!g.edges.get(&s).unwrap().contains(&b));
    
}