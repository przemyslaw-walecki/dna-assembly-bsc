//! Moduł odpowiedzialny za budowę i czyszczenie grafu de Bruijna opartego na k-merach.
//!
//! Graf reprezentuje:
//! - węzły jako k-mery zakodowane w `u128`,
//! - krawędzie jako przejścia (k+1)-merów,
//! - pokrycie węzłów (KC — *k-mer coverage*),
//! - pokrycie krawędzi (EC — *edge coverage*).
//!
//! Udostępnia mechanizmy:
//! - konstrukcji grafu z listy odczytów,
//! - filtrowania niskiego pokrycia,
//! - usuwania krótkich zakończeń (tips),
//! - wykrywania i usuwania prostych bąbli (alternatywnych ścieżek),
//! - eksportu grafu do formatu GFA z informacjami o pokryciu.

use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::collections::{HashSet as StdHashSet, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};

/// Skierowany graf de Bruijna oparty na k-merach.
///
/// Zawiera:
/// - listę krawędzi (następników),
/// - stopnie wejścia/wyjścia,
/// - pokrycie węzłów i krawędzi.
pub struct KmerGraph {
    /// Długość k-mera.
    pub k: usize,
    /// Lista następników dla każdego węzła.
    pub edges: HashMap<u128, HashSet<u128>>,
    /// Liczba krawędzi wchodzących na węzeł.
    pub in_degree: HashMap<u128, usize>,
    /// Liczba krawędzi wychodzących z węzła.
    pub out_degree: HashMap<u128, usize>,
    /// Pokrycie krawędzi (k+1-mery).
    pub edge_counts: HashMap<(u128, u128), usize>,
    /// Pokrycie węzłów (k-mery).
    pub node_coverage: HashMap<u128, usize>,
}

impl KmerGraph {
    /// Tworzy pusty graf o zadanej długości k-mera.
    pub fn new(k: usize) -> Self {
        Self {
            k,
            edges: HashMap::default(),
            in_degree: HashMap::default(),
            out_degree: HashMap::default(),
            edge_counts: HashMap::default(),
            node_coverage: HashMap::default(),
        }
    }

    /// Koduje k-mer (np. "ACGT") do postaci `u128`.
    ///
    /// A=0, C=1, G=2, T=3 — każda zasada kodowana na 2 bitach.
    ///
    /// # Panika
    /// W przypadku napotkania niepoprawnego znaku.
    pub fn encode(&self, s: &str) -> u128 {
        let mut c: u128 = 0;
        for b in s.bytes() {
            let v = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("Niepoprawna zasada: {}", b as char),
            };
            c = (c << 2) | v;
        }
        c
    }

    /// Dekoduje `u128` z powrotem do k-mera.
    pub fn decode(&self, mut code: u128) -> String {
        let mut buf = vec!['A'; self.k];
        for i in (0..self.k).rev() {
            buf[i] = match code & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            code >>= 2;
        }
        buf.into_iter().collect()
    }

    /// Buduje graf z listy odczytów.
    ///
    /// Zlicza:
    /// - pokrycie k-merów,
    /// - pokrycie (k+1)-merów,
    /// - stopnie wejścia/wyjścia.
    pub fn build(&mut self, reads: &[String]) {
        for r in reads {
            if r.len() < self.k {
                continue;
            }

            // Zliczanie k-merów (pokrycie węzłów)
            for i in 0..=r.len() - self.k {
                let u = self.encode(&r[i..i + self.k]);
                *self.node_coverage.entry(u).or_default() += 1;
            }

            // Dodawanie krawędzi i zliczanie (k+1)-merów (pokrycie krawędzi)
            if r.len() >= self.k + 1 {
                for i in 0..=r.len() - self.k - 1 {
                    let u = self.encode(&r[i..i + self.k]);
                    let v = self.encode(&r[i + 1..i + 1 + self.k]);

                    *self.edge_counts.entry((u, v)).or_default() += 1;

                    let succs = self.edges.entry(u).or_default();
                    if succs.insert(v) {
                        *self.out_degree.entry(u).or_default() += 1;
                        *self.in_degree.entry(v).or_default() += 1;
                    }
                }
            }
        }

        self.normalize_nodes();
    }

    /// Usuwa krawędzie o pokryciu mniejszym niż `threshold`.
    pub fn filter_low_coverage(&mut self, threshold: usize) {
        let mut to_remove = Vec::new();

        for (&(u, v), &cnt) in &self.edge_counts {
            if cnt < threshold {
                to_remove.push((u, v));
            }
        }

        for (u, v) in to_remove {
            if let Some(succs) = self.edges.get_mut(&u) {
                if succs.remove(&v) {
                    *self.out_degree.get_mut(&u).unwrap() -= 1;
                    *self.in_degree.get_mut(&v).unwrap() -= 1;
                }
            }
            self.edge_counts.remove(&(u, v));
        }
    }

    /// Usuwa krótkie zakończenia (tips) o długości do `max_depth`.
    ///
    /// Zakończeniem jest węzeł o `out_degree = 0`.  
    /// Jeśli prowadzi do niego łańcuch węzłów 1→1 o długości ≤ `max_depth`,
    /// cały łańcuch zostaje usunięty.
    pub fn remove_dead_ends(&mut self, max_depth: Option<usize>) {
        let max_depth = max_depth.unwrap_or(self.k);

        for depth in 1..=max_depth {
            // Budowa mapy poprzedników
            let mut pred: HashMap<u128, HashSet<u128>> = HashMap::default();
            for (&u, succs) in &self.edges {
                for &v in succs {
                    pred.entry(v).or_default().insert(u);
                }
            }

            // Węzły końcowe (tips)
            let tips: Vec<u128> = self
                .out_degree
                .iter()
                .filter_map(|(&n, &od)| if od == 0 { Some(n) } else { None })
                .collect();

            let mut to_remove = HashSet::default();

            for tip in tips {
                let mut chain = vec![tip];
                let mut curr = tip;

                // Idziemy wstecz dopóki występuje unikalny poprzednik 1→1
                for _ in 0..depth {
                    let preds = pred
                        .get(&curr)
                        .unwrap_or(&HashSet::default())
                        .iter()
                        .filter(|&&p| {
                            self.in_degree.get(&p).copied().unwrap_or(0) == 1
                                && self.out_degree.get(&p).copied().unwrap_or(0) == 1
                        })
                        .copied()
                        .collect::<Vec<_>>();

                    if preds.len() != 1 {
                        break;
                    }

                    curr = preds[0];
                    chain.push(curr);
                }

                if chain.len() - 1 <= depth {
                    to_remove.extend(chain);
                }
            }

            if to_remove.is_empty() {
                break;
            }

            // Usuwanie całych łańcuchów zakończeń
            for &n in &to_remove {
                if let Some(succs) = self.edges.remove(&n) {
                    for v in succs {
                        *self.in_degree.get_mut(&v).unwrap() -= 1;
                        self.edge_counts.remove(&(n, v));
                    }
                }
                if let Some(parents) = pred.get(&n) {
                    for &p in parents {
                        if let Some(succs) = self.edges.get_mut(&p) {
                            if succs.remove(&n) {
                                *self.out_degree.get_mut(&p).unwrap() -= 1;
                                self.edge_counts.remove(&(p, n));
                            }
                        }
                    }
                }
                self.in_degree.remove(&n);
                self.out_degree.remove(&n);
                self.node_coverage.remove(&n);
            }
        }

        self.normalize_nodes();
    }

    /// Usuwa proste bąble (dwie alternatywne ścieżki z tego samego źródła).
    ///
    /// Mechanizm:
    /// - BFS do maksymalnej głębokości `max_depth`,
    /// - jeśli dwa różne pathy dochodzą do tego samego węzła,
    ///   krótszy z nich jest usuwany,
    /// - wewnętrzne węzły krótszej ścieżki trafiają do `to_drop`.
    pub fn remove_bubbles(&mut self, max_depth: usize) {
        // Mapa poprzedników
        let mut pred: HashMap<u128, HashSet<u128>> = HashMap::default();
        for (&u, succs) in &self.edges {
            for &v in succs {
                pred.entry(v).or_default().insert(u);
            }
        }

        let mut to_drop: HashSet<u128> = HashSet::default();
        let mut seen_pairs: HashSet<(Vec<u128>, Vec<u128>)> = HashSet::default();

        for (&src, succs) in &self.edges {
            if succs.len() < 2 {
                continue;
            }

            let mut seen: HashMap<u128, Vec<u128>> = HashMap::default();
            let mut q: VecDeque<Vec<u128>> = VecDeque::new();
            q.push_back(vec![src]);

            while let Some(path) = q.pop_front() {
                if path.len() > max_depth {
                    continue;
                }

                let n = *path.last().unwrap();

                if let Some(neighbors) = self.edges.get(&n) {
                    for &nxt in neighbors {
                        if path.contains(&nxt) {
                            continue;
                        }

                        let mut newp = path.clone();
                        newp.push(nxt);

                        if let Some(p2) = seen.get(&nxt) {
                            let key = (newp.clone(), p2.clone());

                            if !seen_pairs.contains(&key) {
                                seen_pairs.insert(key.clone());
                                seen_pairs.insert((p2.clone(), newp.clone()));

                                // Usuwamy węzły krótszej ścieżki
                                let drop_path =
                                    if newp.len() < p2.len() { &newp } else { p2 };

                                to_drop.extend(drop_path.iter().copied());
                            }
                        } else {
                            seen.insert(nxt, newp.clone());
                            q.push_back(newp);
                        }
                    }
                }
            }
        }

        // Usuwanie węzłów należących do krótszych alternatywnych ścieżek
        for &n in &to_drop {
            if let Some(succs) = self.edges.remove(&n) {
                for v in succs {
                    *self.in_degree.get_mut(&v).unwrap() -= 1;
                    self.edge_counts.remove(&(n, v));
                }
            }

            if let Some(parents) = pred.get(&n) {
                for &p in parents {
                    if let Some(succs) = self.edges.get_mut(&p) {
                        if succs.remove(&n) {
                            *self.out_degree.get_mut(&p).unwrap() -= 1;
                            self.edge_counts.remove(&(p, n));
                        }
                    }
                }
            }

            self.in_degree.remove(&n);
            self.out_degree.remove(&n);
            self.node_coverage.remove(&n);
        }

        self.normalize_nodes();
    }

    /// Zwraca całkowitą liczbę krawędzi w grafie.
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|s| s.len()).sum()
    }

    /// Upewnia się, że każdy węzeł posiada wpisy KC, stopni wejścia i wyjścia.
    ///
    /// Jest to wymagane po operacjach czyszczących, które mogły usunąć niektóre dane.
    fn normalize_nodes(&mut self) {
        let all_nodes: StdHashSet<u128> = self
            .edges
            .keys()
            .copied()
            .chain(self.in_degree.keys().copied())
            .chain(self.out_degree.keys().copied())
            .chain(self.node_coverage.keys().copied())
            .collect();

        for n in all_nodes {
            self.edges.entry(n).or_default();
            self.in_degree.entry(n).or_default();
            self.out_degree.entry(n).or_default();
            self.node_coverage.entry(n).or_default();
        }
    }

    /// Zapisuje graf do pliku w formacie GFA.
    ///
    /// Dopisuje:
    /// - `KC:i:<wartość>` — pokrycie węzła,
    /// - `EC:i:<wartość>` — pokrycie krawędzi.
    pub fn write_gfa(&self, path: &str) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "H\tVN:Z:1.0")?;

        // Zbieramy wszystkie węzły, aby nie pominąć izolowanych k-merów.
        let mut all_nodes = StdHashSet::new();

        for (&u, succs) in &self.edges {
            all_nodes.insert(u);
            for &v in succs {
                all_nodes.insert(v);
            }
        }

        for &n in self.in_degree.keys() {
            all_nodes.insert(n);
        }
        for &n in self.out_degree.keys() {
            all_nodes.insert(n);
        }
        for &n in self.node_coverage.keys() {
            all_nodes.insert(n);
        }

        let mut nodes: Vec<_> = all_nodes.into_iter().collect();
        nodes.sort_unstable();

        // Sekcje węzłów S
        for &n in &nodes {
            let seq = self.decode(n);
            let kc = self.node_coverage.get(&n).copied().unwrap_or(0);
            writeln!(w, "S\t{}\t{}\tKC:i:{}", n, seq, kc)?;
        }

        // Sekcje krawędzi L
        for (&u, succs) in &self.edges {
            for &v in succs {
                let ec = self.edge_counts.get(&(u, v)).copied().unwrap_or(0);
                writeln!(w, "L\t{}\t+\t{}\t+\t0M\tEC:i:{}", u, v, ec)?;
            }
        }

        Ok(())
    }
}
