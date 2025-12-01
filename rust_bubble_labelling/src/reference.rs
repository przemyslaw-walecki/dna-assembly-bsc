//! # Moduł referencyjny
//!
//! Zawiera funkcje pomocnicze związane z sekwencją referencyjną:
//! - reverse complement,
//! - ładowanie FASTA,
//! - wyszukiwanie pozycji k-merów,
//! - minimal unique extension (MUE).

use std::fs::File;
use std::io::{BufRead, BufReader};

#[inline]
pub fn revcomp(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.bytes().rev() {
        out.push(match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        });
    }
    out
}

/// Pojedynczy chromosom.
#[derive(Clone)]
pub struct Chrom {
    pub name: String,
    pub seq: String,
    pub rc_seq: String,
}

/// Struktura referencji, obsługująca wiele chromosomów.
pub struct RefText {
    pub chroms: Vec<Chrom>,
}

impl RefText {
    /// Ładowanie pliku FASTA do pamięci.
    pub fn load_fasta(path: &str) -> Self {
        let f = File::open(path).expect("Nie można otworzyć FASTA");
        let rdr = BufReader::new(f);

        let mut name = String::new();
        let mut seq = String::new();
        let mut chroms_raw = Vec::<(String, String)>::new();

        for line in rdr.lines().flatten() {
            if line.starts_with('>') {
                if !seq.is_empty() {
                    chroms_raw.push((name.clone(), seq.clone()));
                    seq.clear();
                }
                name = line[1..].trim().to_string();
            } else {
                seq.push_str(line.trim());
            }
        }
        if !seq.is_empty() {
            chroms_raw.push((name, seq));
        }

        let chroms = chroms_raw
            .into_iter()
            .map(|(n, s)| Chrom {
                name: n,
                rc_seq: revcomp(&s),
                seq: s,
            })
            .collect();

        Self { chroms }
    }
}

/// Struktura pozycji k-mera na referencji.
#[derive(Clone, Copy)]
pub struct Pos {
    pub chrom: usize,
    pub pos: usize,
    pub rev: bool,
}

/// Wyszukiwanie wszystkich wystąpień k-mera (forward + RC).
pub fn all_kmer_pos(reftext: &RefText, kmer: &str) -> Vec<Pos> {
    fn scan(hay: &str, needle: &str, chrom: usize, rev: bool, out: &mut Vec<Pos>) {
        if needle.is_empty() || hay.len() < needle.len() {
            return;
        }
        let nb = needle.as_bytes();
        let hb = hay.as_bytes();
        let m = nb.len();

        let mut i = 0usize;
        while i + m <= hb.len() {
            if &hb[i..i + m] == nb {
                out.push(Pos { chrom, pos: i, rev });
            }
            i += 1;
        }
    }

    let mut out = Vec::new();
    let rc = revcomp(kmer);

    for (i, c) in reftext.chroms.iter().enumerate() {
        scan(&c.seq, kmer, i, false, &mut out);
        scan(&c.rc_seq, &rc, i, true, &mut out);
    }

    out
}

/// Minimal Unique Extension (MUE) — rozszerza k-mer, aż stanie się unikatowy.
pub fn minimal_unique_extension(
    reftext: &RefText,
    base_kmer: &str,
    extend_right: bool,
    max_ext: usize,
) -> Option<(usize, Pos)> {
    use std::collections::HashMap;

    let mut occ = all_kmer_pos(reftext, base_kmer);
    if occ.is_empty() {
        return None;
    }
    if occ.len() == 1 {
        return Some((0, occ[0]));
    }

    for ext in 1..=max_ext {
        let mut groups: HashMap<(usize, u8, bool), Vec<Pos>> = HashMap::new();

        for p in occ.iter().copied() {
            let c = &reftext.chroms[p.chrom];

            let next_base = if !p.rev {
                if extend_right {
                    c.seq.as_bytes().get(p.pos + base_kmer.len() + (ext - 1))
                } else {
                    p.pos.checked_sub(ext)
                        .and_then(|i| c.seq.as_bytes().get(i))
                }
            } else {
                if extend_right {
                    c.rc_seq.as_bytes().get(p.pos + base_kmer.len() + (ext - 1))
                } else {
                    p.pos.checked_sub(ext)
                        .and_then(|i| c.rc_seq.as_bytes().get(i))
                }
            };

            if let Some(&b) = next_base {
                groups.entry((p.chrom, b, p.rev))
                    .or_default()
                    .push(p);
            }
        }

        let mut uniq = Vec::<Pos>::new();
        for v in groups.values() {
            if v.len() == 1 {
                uniq.push(v[0]);
            }
        }

        if uniq.len() == 1 {
            return Some((ext, uniq[0]));
        }

        occ = groups.values().flat_map(|v| v.clone()).collect();
        if occ.len() == 1 {
            return Some((ext, occ[0]));
        }
    }

    None
}
