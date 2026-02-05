use std::path::PathBuf;

use rust_basic_dna_assembler::{Assembler, FastqParser, KmerGraph};

fn fixture(path: &str) -> String {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/fastq");
    p.push(path);
    p.to_string_lossy().to_string()
}

fn load_reads_from_fixture(path: &str) -> Vec<String> {
    let parser = FastqParser::new(fixture(path));
    parser.load_reads().expect("failed to load fixture FASTQ")
}

fn build_graph_from_two_fastqs(k: usize, f1: &str, f2: &str) -> KmerGraph {
    let mut reads = Vec::new();
    reads.extend(load_reads_from_fixture(f1));
    reads.extend(load_reads_from_fixture(f2));

    let mut g = KmerGraph::new(k);
    g.build(&reads);
    g
}

fn rotations(s: &str) -> Vec<String> {
    let bytes = s.as_bytes();
    let n = bytes.len();
    (0..n)
        .map(|i| {
            let mut out = Vec::with_capacity(n);
            out.extend_from_slice(&bytes[i..]);
            out.extend_from_slice(&bytes[..i]);
            String::from_utf8(out).unwrap()
        })
        .collect()
}

#[test]
fn pipeline_linear_fixture_produces_single_unitig() {
    // linear.fastq: ACGT
    let g = build_graph_from_two_fastqs(3, "linear.fastq", "linear.fastq");
    let asm = Assembler::new(&g);

    let unitigs = asm.assemble_unitigs();

    assert_eq!(unitigs.len(), 1);
    assert_eq!(unitigs[0], "ACGT");
}

#[test]
fn pipeline_branching_fixture_produces_two_unitigs() {
    // branching.fastq: AAATG i AAACG
    let g = build_graph_from_two_fastqs(3, "branching.fastq", "branching.fastq");
    let asm = Assembler::new(&g);

    let mut unitigs = asm.assemble_unitigs();
    unitigs.sort();

    assert_eq!(unitigs.len(), 2);
    assert!(unitigs.contains(&"AAACG".to_string()));
    assert!(unitigs.contains(&"AAATG".to_string()));
}

#[test]
fn pipeline_cycle_fixture_produces_single_cycle_unitig() {
    // cycle.fastq: ATGATG -> czysty cykl 1->1 na k=3
    let g = build_graph_from_two_fastqs(3, "cycle.fastq", "cycle.fastq");
    let asm = Assembler::new(&g);

    let unitigs = asm.assemble_unitigs();
    assert_eq!(unitigs.len(), 1);

    let got = &unitigs[0];
    // dla k=3 i cyklu 3-węzłowego unitig ma długość 5 (ATGAT) z rotacją
    assert_eq!(got.len(), 5);

    let expected = "ATGAT";
    let rots = rotations(expected);
    assert!(rots.contains(got), "unexpected cycle unitig: {got}");
}
