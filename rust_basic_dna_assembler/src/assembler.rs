//! Moduł odpowiedzialny za składanie kontigów w assemblerze DNA.
//!
//! Udostępnia strukturę [`Assembler`], która generuje **unitigi** (maksymalne
//! ścieżki 1->1), w tym cykle o stopniu wejścia i wyjścia równym 1.  
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
/// strukturze [`KmerGraph`]. Przechodzi ścieżki 1->1 i składa odpowiadające im
/// sekwencje DNA. Obsługiwane są również cykle czysto 1->1.
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

    /// Składa kontigi metodą opartą na klasycznych ścieżkach 1->1.
    ///
    /// Metoda utrzymana jako nieużywana: zachowana w celu porównawczym i
    /// dokumentacyjnym.
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

    /// Składa wszystkie unitigi, w tym cykle o stopniach 1->1.
    ///
    /// Działanie obejmuje:
    /// - znajdowanie ścieżek rozpoczynających się w węzłach nie-1->1,
    /// - rozszerzanie ich dopóki trwa sekwencja 1->1,
    /// - osobną obsługę cykli, które nie mają początku ani końca,
    /// - unikanie powtórnego odwiedzania tych samych krawędzi.
    ///
    /// Zwraca:
    /// - wektor sekwencji DNA reprezentujących unitigi.
    pub fn assemble_unitigs(&self) -> Vec<String> {
        let mut unitigs = Vec::new();
        let mut visited: HashSet<(u128, u128)> = HashSet::new();

        // --- Krok 1: ścieżki zaczynające się w węzłach nie-1->1 ----------------
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

        // --- Krok 2: obsługa czystych cykli 1->1 ------------------------------
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer_graph::KmerGraph;

    fn build_graph(k: usize, reads: &[&str]) -> KmerGraph {
        let mut g = KmerGraph::new(k);
        let v = reads.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        g.build(&v);
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
    fn assemble_unitigs_linear_path_single_contig() {
        // k=3, read: ACGT => ACG -> CGT
        let g = build_graph(3, &["ACGT"]);
        let asm = Assembler::new(&g);

        let unitigs = asm.assemble_unitigs();
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0], "ACGT");
    }

    #[test]
    fn assemble_contigs_linear_path_single_contig() {
        let g = build_graph(3, &["ACGT"]);
        let asm = Assembler::new(&g);

        let contigs = asm.assemble_contigs();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0], "ACGT");
    }

    #[test]
    fn assemble_unitigs_branching_produces_two_unitigs() {
        // Graf:
        // AAA -> AAT -> ATG  => "AAATG"
        // AAA -> AAC -> ACG  => "AAACG"
        //
        // Osiągamy to przez dwa odczyty długości 5 (k=3).
        let g = build_graph(3, &["AAATG", "AAACG"]);
        let asm = Assembler::new(&g);

        let mut unitigs = asm.assemble_unitigs();
        unitigs.sort();

        assert_eq!(unitigs.len(), 2);
        assert!(unitigs.contains(&"AAATG".to_string()));
        assert!(unitigs.contains(&"AAACG".to_string()));
    }

    #[test]
    fn assemble_unitigs_pure_1to1_cycle_detected() {
        // Czysty cykl 1->1 na k=3:
        // ATG -> TGA -> GAT -> ATG
        //
        // Odczyt "ATGATG" daje kmery: ATG, TGA, GAT, ATG i domyka cykl.
        let g = build_graph(3, &["ATGATG"]);
        let asm = Assembler::new(&g);

        let unitigs = asm.assemble_unitigs();

        // Dla czystego cyklu powinniśmy dostać dokładnie jeden unitig.
        assert_eq!(unitigs.len(), 1);

        // Reprezentacja cyklu może zacząć się w dowolnym miejscu (rotacja),
        // ale długość powinna wynosić k + (liczba_węzłów - 1) = 3 + 2 = 5.
        let u = &unitigs[0];
        assert_eq!(u.len(), 5);

        // Jedna z rotacji "ATGAT" (zamknięty cykl bez duplikowania ostatniej zasady)
        // powinna pasować.
        let allowed = ["ATGAT", "TGATG", "GATGA"];
        assert!(
            allowed.contains(&u.as_str()),
            "unexpected cycle unitig: {u}"
        );
    }

    #[test]
    fn reconstruct_respects_cycle_flag() {
        // W tym teście celowo używamy prywatnej metody `reconstruct`,
        // co jest dozwolone, bo testy są w module potomnym.
        let g = build_graph(3, &["ATGATG"]);
        let asm = Assembler::new(&g);

        let atg = g.encode("ATG");
        let tga = g.encode("TGA");
        let gat = g.encode("GAT");

        // Ścieżka bez domknięcia cyklu: ATG -> TGA -> GAT
        let linear = asm.reconstruct(&[atg, tga, gat], false);
        assert_eq!(linear, "ATGAT"); // 3 + 2

        // Ścieżka z domknięciem cyklu: ATG -> TGA -> GAT -> ATG
        // Przy is_cycle=true nie doklejamy ostatniej zasady dla powtórzonego startu,
        // więc wynik powinien być taki sam jak wyżej.
        let cycle = asm.reconstruct(&[atg, tga, gat, atg], true);
        assert_eq!(cycle, "ATGAT");
    }

    #[test]
    fn assemble_unitigs_does_not_duplicate_edges_across_outputs() {
        // Weryfikacja, że visited na krawędziach działa: dla liniowego grafu nie powstaną duplikaty.
        let g = build_graph(3, &["ACGT", "ACGT"]); // ten sam read 2x zwiększa coverage, ale nie powinien tworzyć 2 unitigów
        let asm = Assembler::new(&g);

        let unitigs = asm.assemble_unitigs();
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0], "ACGT");
    }
}
