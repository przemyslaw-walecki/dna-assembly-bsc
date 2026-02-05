//! # Moduł GFA
//!
//! Ten moduł odpowiada za:
//! - wczytywanie grafu z pliku GFA,
//! - przechowywanie segmentów (`S`) oraz krawędzi (`L`),
//! - budowanie uproszczonych struktur następników i poprzedników.
//!
//! Używany przez pipeline etykietowania bąbli hiperstrukturalnych.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Segment (`S` w GFA) – węzeł grafu.
#[derive(Debug)]
pub struct Segment {
    pub id: u128,
    pub seq: String,
    pub cov: usize, // KC
}

/// Krawędź (`L` w GFA).
#[derive(Debug, Clone)]
pub struct Link {
    pub from: u128,
    pub to: u128,
    pub cov: usize, // EC
}

/// Wczytuje graf z pliku GFA oraz zwraca:
/// - mapę segmentów,
/// - listę krawędzi,
/// - słownik pokryć EC.
pub fn parse_gfa(
    path: &str,
) -> (
    HashMap<u128, Segment>,
    Vec<Link>,
    HashMap<(u128, u128), usize>,
) {
    let f = File::open(path).expect("Nie można otworzyć pliku GFA");
    let rdr = BufReader::new(f);

    let mut segs = HashMap::<u128, Segment>::new();
    let mut links = Vec::<Link>::new();
    let mut edge_cov = HashMap::<(u128, u128), usize>::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();

        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 {
                    continue;
                }
                let id: u128 = cols[1]
                    .parse()
                    .unwrap_or_else(|_| panic!("Niepoprawny numeric id segmentu: {}", cols[1]));
                let seq = cols[2].to_string();
                let cov = cols
                    .iter()
                    .find_map(|f| {
                        let lf = f.to_ascii_lowercase();
                        lf.strip_prefix("kc:i:").and_then(|x| x.parse::<usize>().ok())
                    })
                    .unwrap_or(0);
                segs.insert(id, Segment { id, seq, cov });
            }
            Some("L") => {
                if cols.len() < 6 {
                    continue;
                }
                let from: u128 = cols[1]
                    .parse()
                    .unwrap_or_else(|_| panic!("Niepoprawne id 'from': {}", cols[1]));
                let to: u128 = cols[3]
                    .parse()
                    .unwrap_or_else(|_| panic!("Niepoprawne id 'to': {}", cols[3]));

                let cov = cols
                    .iter()
                    .find_map(|f| {
                        let lf = f.to_ascii_lowercase();
                        lf.strip_prefix("ec:i:").and_then(|x| x.parse::<usize>().ok())
                    })
                    .unwrap_or(0);

                links.push(Link { from, to, cov });
                edge_cov.insert((from, to), cov);
            }
            _ => {}
        }
    }

    (segs, links, edge_cov)
}

/// Buduje mapy następników (`succ`) i poprzedników (`pred`) z listy krawędzi.
pub fn build_adj(
    segs: &HashMap<u128, Segment>,
    links: &[Link],
) -> (
    HashMap<u128, Vec<u128>>,
    HashMap<u128, Vec<u128>>,
) {
    let mut succ: HashMap<u128, Vec<u128>> =
        segs.keys().map(|&k| (k, Vec::new())).collect();
    let mut pred: HashMap<u128, Vec<u128>> =
        segs.keys().map(|&k| (k, Vec::new())).collect();

    for l in links {
        if let Some(v) = succ.get_mut(&l.from) {
            v.push(l.to);
        }
        if let Some(v) = pred.get_mut(&l.to) {
            v.push(l.from);
        }
    }

    (succ, pred)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        // dość unikatowo bez dodatkowych crate’ów
        let pid = std::process::id();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("{}_{}_{}_{}.tmp", name, pid, nanos, "gfa"));
        p
    }

    #[test]
    fn parse_gfa_reads_segments_links_and_coverages_case_insensitive() {
        let p = tmp_path("parse_gfa_basic");

        // Uwaga: parser bierze:
        // S: cols[1]=id, cols[2]=seq, pole KC:I:... (case-insensitive)
        // L: cols[1]=from, cols[3]=to, pole EC:i:... (case-insensitive)
        let gfa = "\
S\t1\tACGT\tKC:i:7
S\t2\tCGTA\tkc:I:9
S\t3\tGTAA
L\t1\t+\t2\t+\t0M\tEC:I:5
L\t2\t+\t3\t+\t0M\tec:i:11
";
        fs::write(&p, gfa).unwrap();

        let (segs, links, edge_cov) = parse_gfa(p.to_str().unwrap());

        assert_eq!(segs.len(), 3);
        assert_eq!(segs[&1].seq, "ACGT");
        assert_eq!(segs[&1].cov, 7);
        assert_eq!(segs[&2].cov, 9);
        assert_eq!(segs[&3].cov, 0); // brak KC => 0

        assert_eq!(links.len(), 2);
        assert_eq!(links[0].from, 1);
        assert_eq!(links[0].to, 2);
        assert_eq!(links[0].cov, 5);
        assert_eq!(links[1].cov, 11);

        assert_eq!(edge_cov[&(1, 2)], 5);
        assert_eq!(edge_cov[&(2, 3)], 11);

        let _ = fs::remove_file(&p);
    }

    #[test]
    fn build_adj_creates_succ_and_pred_for_all_nodes() {
        let mut segs = HashMap::<u128, Segment>::new();
        segs.insert(1, Segment { id: 1, seq: "AAA".into(), cov: 1 });
        segs.insert(2, Segment { id: 2, seq: "AAT".into(), cov: 2 });
        segs.insert(3, Segment { id: 3, seq: "ATG".into(), cov: 3 });

        let links = vec![
            Link { from: 1, to: 2, cov: 10 },
            Link { from: 2, to: 3, cov: 20 },
        ];

        let (succ, pred) = build_adj(&segs, &links);

        // wszystkie węzły obecne jako klucze
        assert!(succ.contains_key(&1) && succ.contains_key(&2) && succ.contains_key(&3));
        assert!(pred.contains_key(&1) && pred.contains_key(&2) && pred.contains_key(&3));

        assert_eq!(succ[&1], vec![2]);
        assert_eq!(succ[&2], vec![3]);
        assert!(succ[&3].is_empty());

        assert!(pred[&1].is_empty());
        assert_eq!(pred[&2], vec![1]);
        assert_eq!(pred[&3], vec![2]);
    }

    #[test]
    #[should_panic]
    fn parse_gfa_panics_on_non_numeric_segment_id() {
        let p = tmp_path("parse_gfa_bad_id");
        let gfa = "S\tabc\tACGT\tKC:i:1\n";
        fs::write(&p, gfa).unwrap();

        let _ = parse_gfa(p.to_str().unwrap());
    }
}
