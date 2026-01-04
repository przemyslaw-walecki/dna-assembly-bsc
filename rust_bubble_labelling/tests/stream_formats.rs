mod common;
use common::{make_tmp_dir, write_file, run_bin, read_lines};

#[test]
fn reads_multiple_json_values_and_multiple_shapes() {
    let dir = make_tmp_dir("stream_formats");

    // Segmenty muszą istnieć, bo start/end są używane do seqs itd.
    let gfa_txt = "\
S\t1\tACGT\tKC:i:10
S\t2\tCGTA\tKC:i:12
S\t3\tGTAA\tKC:i:8
";
    let gfa = write_file(&dir, "graph.gfa", gfa_txt);

    // Dwa JSON-y jeden po drugim:
    // 1) tablica Bubble
    // 2) obiekt mapujący do BubbleChain
    let bubbles_stream = r#"
[{"id":1,"ends":["1","2"],"inside":[]}]
{"x":{"bubbles":[{"id":2,"ends":["2","3"],"inside":[]}]}}
"#;
    let bubbles = write_file(&dir, "bubbles.jsonl", bubbles_stream);

    let reference = write_file(&dir, "ref.fa", ">chr1\nACGTACGTACGTACGT\n");
    let out = dir.join("out.jsonl");

    let (code, _stdout, _stderr) = run_bin(&[
        gfa.to_str().unwrap(),
        bubbles.to_str().unwrap(),
        reference.to_str().unwrap(),
        out.to_str().unwrap(),
    ]);

    assert_eq!(code, 0);

    let lines = read_lines(&out);
    assert_eq!(lines.len(), 2);

    let v1: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
    let v2: serde_json::Value = serde_json::from_str(&lines[1]).unwrap();

    assert_eq!(v1["bubble_id"].as_u64().unwrap(), 1);
    assert_eq!(v2["bubble_id"].as_u64().unwrap(), 2);

    // Dla takiego grafu bez linków nadal spodziewamy się no_paths
    assert_eq!(v1["label_reason"].as_str().unwrap(), "no_paths");
    assert_eq!(v2["label_reason"].as_str().unwrap(), "no_paths");
}
