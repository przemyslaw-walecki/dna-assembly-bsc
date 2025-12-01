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
