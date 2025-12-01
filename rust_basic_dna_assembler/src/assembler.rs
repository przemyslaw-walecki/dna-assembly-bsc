//! Moduł odpowiedzialny za składanie kontigów w assemblerze DNA.
//!
//! Udostępnia strukturę [`Assembler`], która generuje **unitigi** (maksymalne
//! ścieżki 1→1), w tym cykle o stopniu wejścia i wyjścia równym 1.  
//! Implementacja obejmuje dwa algorytmy:
//! - `assemble_contigs()` — metoda starsza, zachowana jedynie jako kod nieużywany,
//! - `assemble_unitigs()` — domyślne tworzenie unitigów z obsługą cykli.
//!
//! Rekonstrukcja sekwencji opiera się na dekodowaniu k-merów przechowywanych
//! w strukturze [`KmerGraph`].

use crate::kmer_graph::KmerGraph;
use std::collections::HashSet;

/// Assembler generujący kontigi (unitigi) na podstawie grafu de Bruijna.
///
/// Assembler nie modyfikuje grafu – działa tylko na gotowej, oczyszczonej
/// strukturze [`KmerGraph`]. Przechodzi ścieżki 1→1 i składa odpowiadające im
/// sekwencje DNA. Obsługiwane są również cykle czysto 1→1.
pub struct Assembler<'g> {
    /// Współdzielona referencja do grafu k-merów.
    pub graph: &'g KmerGraph,
}

impl<'g> Assembler<'g> {
    /// Tworzy nowy assembler korzystający z przekazanego grafu.
    pub fn new(graph: &'g KmerGraph) -> Self {
        Self { graph }
    }

    // -------------------------------------------------------------------------
    //  STARA METODA – pozostawiona jako kod martwy
    // -------------------------------------------------------------------------

    /// Składa kontigi metodą opartą na klasycznych ścieżkach 1→1.
    ///
    /// Metoda utrzymana jako nieużywana: zachowana w celu porównawczym i
    /// dokumentacyjnym. Domyślnym algorytmem produkcji kontigów jest
    /// [`assemble_unitigs`].
    pub fn assemble_contigs(&self) -> Vec<String> {
        let mut contigs = Vec::new();
        let mut visited = HashSet::new();

        for (&u, succs) in &self.graph.edges {
            let indeg = *self.graph.in_degree.get(&u).unwrap_or(&0);
            let outdeg = *self.graph.out_degree.get(&u).unwrap_or(&0);

            if indeg != 1 || outdeg != 1 {
                for &v in succs {
                    if !visited.insert((u, v)) {
                        continue;
                    }

                    let mut path = vec![u, v];
                    let mut curr = v;

                    while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                        if let Some(nxt) = self.first_succ(curr) {
                            if !visited.insert((curr, nxt)) {
                                break;
                            }
                            curr = nxt;
                            path.push(curr);
                        } else {
                            break;
                        }
                    }

                    contigs.push(self.reconstruct(&path, false));
                }
            }
        }
        contigs
    }

    /// Składa wszystkie unitigi, w tym cykle o stopniach 1→1.
    ///
    /// Działanie obejmuje:
    /// - znajdowanie ścieżek rozpoczynających się w węzłach nie-1→1,
    /// - rozszerzanie ich dopóki trwa sekwencja 1→1,
    /// - osobną obsługę cykli, które nie mają początku ani końca,
    /// - unikanie powtórnego odwiedzania tych samych krawędzi.
    ///
    /// Zwraca:
    /// - wektor sekwencji DNA reprezentujących unitigi.
    pub fn assemble_unitigs(&self) -> Vec<String> {
        let mut unitigs = Vec::new();
        let mut visited: HashSet<(u128, u128)> = HashSet::new();

        // --- Krok 1: ścieżki zaczynające się w węzłach nie-1→1 ----------------
        for (&u, succs) in &self.graph.edges {
            let indeg = self.deg_in(u);
            let outdeg = self.deg_out(u);

            if indeg != 1 || outdeg != 1 {
                for &v in succs {
                    if !visited.insert((u, v)) {
                        continue;
                    }

                    let mut path = vec![u, v];
                    let mut curr = v;

                    while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                        if let Some(nxt) = self.first_succ(curr) {
                            if !visited.insert((curr, nxt)) {
                                break;
                            }
                            curr = nxt;
                            path.push(curr);
                        } else {
                            break;
                        }
                    }

                    unitigs.push(self.reconstruct(&path, false));
                }
            }
        }

        // --- Krok 2: obsługa czystych cykli 1→1 ------------------------------
        for (&u, succs) in &self.graph.edges {
            if self.deg_in(u) != 1 || self.deg_out(u) != 1 {
                continue;
            }

            for &v in succs {
                if visited.contains(&(u, v)) {
                    continue;
                }

                let mut path = vec![u, v];
                visited.insert((u, v));
                let mut curr = v;

                loop {
                    if self.deg_in(curr) != 1 || self.deg_out(curr) != 1 {
                        while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                            if let Some(nxt) = self.first_succ(curr) {
                                if !visited.insert((curr, nxt)) {
                                    break;
                                }
                                curr = nxt;
                                path.push(curr);
                            } else {
                                break;
                            }
                        }
                        unitigs.push(self.reconstruct(&path, false));
                        break;
                    }

                    let Some(nxt) = self.first_succ(curr) else {
                        unitigs.push(self.reconstruct(&path, false));
                        break;
                    };

                    if curr == u && nxt == v {
                        path.push(u);
                        unitigs.push(self.reconstruct(&path, true));
                        break;
                    }

                    if !visited.insert((curr, nxt)) {
                        break;
                    }

                    curr = nxt;
                    path.push(curr);
                }
            }
        }

        unitigs
    }

    // -------------------------------------------------------------------------
    //  METODY POMOCNICZE
    // -------------------------------------------------------------------------

    #[inline]
    fn deg_in(&self, u: u128) -> usize {
        *self.graph.in_degree.get(&u).unwrap_or(&0)
    }

    #[inline]
    fn deg_out(&self, u: u128) -> usize {
        *self.graph.out_degree.get(&u).unwrap_or(&0)
    }

    #[inline]
    fn first_succ(&self, u: u128) -> Option<u128> {
        self.graph.edges.get(&u).and_then(|s| s.iter().next().copied())
    }

    /// Rekonstruuje sekwencję DNA na podstawie ścieżki k-merów.
    ///
    /// Jeśli `is_cycle = true`, ostatni element ścieżki jest powtórzeniem
    /// węzła początkowego – nie należy doklejać jego ostatniej zasady.
    fn reconstruct(&self, path: &[u128], is_cycle: bool) -> String {
        if path.is_empty() {
            return String::new();
        }
        let mut seq = self.graph.decode(path[0]);
        let end = if is_cycle { path.len() - 1 } else { path.len() };
        for &code in &path[1..end] {
            let kmer = self.graph.decode(code);
            seq.push_str(&kmer[kmer.len() - 1..]);
        }
        seq
    }
}
