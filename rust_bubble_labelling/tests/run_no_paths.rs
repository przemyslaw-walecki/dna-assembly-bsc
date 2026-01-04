mod common;
use common::{make_tmp_dir, write_file, run_bin, read_lines};

#[test]
fn generates_record_with_no_paths_reason_when_disconnected() {
    let dir = make_tmp_dir("no_paths");

    // 2 segmenty, brak linków => brak ścieżek między 1 a 2
    let gfa_txt = "\
S\t1\tACGT\tKC:i:10
S\t2\tCGTA\tKC:i:12
";
    let bubbles_txt = r#"[{"id":1,"ends":["1","2"],"inside":[]}]"#;
    let ref_txt = ">chr1\nACGTACGTACGT\n";

    let gfa = write_file(&dir, "graph.gfa", gfa_txt);
    let bubbles = write_file(&dir, "bubbles.json", bubbles_txt);
    let reference = write_file(&dir, "ref.fa", ref_txt);
    let out = dir.join("out.jsonl");

    let (code, _stdout, _stderr) = run_bin(&[
        gfa.to_str().unwrap(),
        bubbles.to_str().unwrap(),
        reference.to_str().unwrap(),
        out.to_str().unwrap(),
    ]);

    assert_eq!(code, 0);

    let lines = read_lines(&out);
    assert_eq!(lines.len(), 1);

    let v: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
    assert_eq!(v["bubble_id"].as_u64().unwrap(), 1);
    assert_eq!(v["k"].as_u64().unwrap(), 4);
    assert_eq!(v["start_seq"].as_str().unwrap(), "ACGT");
    assert_eq!(v["end_seq"].as_str().unwrap(), "CGTA");
    assert_eq!(v["label_reason"].as_str().unwrap(), "no_paths");
    assert!(v["label_path"].is_null());
}
