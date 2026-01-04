mod common;
use common::{make_tmp_dir, write_file, run_bin};

#[test]
fn empty_gfa_exits_nonzero() {
    let dir = make_tmp_dir("empty_graph");

    let gfa = write_file(&dir, "graph.gfa", ""); // pusty
    let bubbles = write_file(&dir, "bubbles.json", r#"[{"id":1,"ends":["1","2"],"inside":[]}]"#);
    let reference = write_file(&dir, "ref.fa", ">chr1\nACGTACGT\n");
    let out = dir.join("out.jsonl");

    let (code, _stdout, stderr) = run_bin(&[
        gfa.to_str().unwrap(),
        bubbles.to_str().unwrap(),
        reference.to_str().unwrap(),
        out.to_str().unwrap(),
    ]);

    assert_ne!(code, 0);
    assert!(stderr.contains("Graf nie zawiera segmentów"));
}
