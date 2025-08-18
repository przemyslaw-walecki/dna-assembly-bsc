use serde::Deserialize;
use serde_json::Value;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

// ====== Small utils ======

fn revcomp(s: &str) -> String {
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
fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>().unwrap_or_else(|_| panic!("Non-numeric id: {s}"))
}

// ====== Low-k graph & bubble handling (candidate enumeration) ======

#[derive(Debug)]
struct Segment { id: u128, seq: String }
#[derive(Debug, Clone)]
struct Link { from: u128, to: u128 }

fn parse_gfa(path: &str) -> (HashMap<u128, Segment>, Vec<Link>) {
    let f = File::open(path).expect("open GFA");
    let rdr = BufReader::new(f);
    let mut segs = HashMap::new();
    let mut links = Vec::new();
    for line in rdr.lines().flatten() {
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 { continue; }
                let id: u128 = cols[1].parse().unwrap();
                let seq = cols[2].to_string();
                segs.insert(id, Segment { id, seq });
            }
            Some("L") => {
                if cols.len() < 6 { continue; }
                let from: u128 = cols[1].parse().unwrap();
                let to: u128 = cols[3].parse().unwrap();
                links.push(Link { from, to });
            }
            _ => {}
        }
    }
    (segs, links)
}
fn build_adj(
    segs: &HashMap<u128, Segment>,
    links: &[Link],
) -> (HashMap<u128, Vec<u128>>, HashMap<u128, Vec<u128>>) {
    let mut succ: HashMap<u128, Vec<u128>> = segs.keys().copied().map(|k| (k, Vec::new())).collect();
    let mut pred: HashMap<u128, Vec<u128>> = segs.keys().copied().map(|k| (k, Vec::new())).collect();
    for l in links {
        if let Some(v) = succ.get_mut(&l.from) { v.push(l.to); }
        if let Some(v) = pred.get_mut(&l.to)   { v.push(l.from); }
    }
    (succ, pred)
}

// compact chains between non-1/1-degree nodes to speed path enumeration
#[derive(Clone)]
struct SuperEdge { from: u128, to: u128, chain: Vec<u128> }
fn compact_with_chains(
    nodes: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
) -> (Vec<SuperEdge>, Vec<u128>) {
    let mut indeg: HashMap<u128, usize> = HashMap::new();
    let mut outdeg: HashMap<u128, usize> = HashMap::new();
    for &n in nodes {
        indeg.insert(n, pred.get(&n).map(|v| v.iter().filter(|x| nodes.contains(x)).count()).unwrap_or(0));
        outdeg.insert(n, succ.get(&n).map(|v| v.iter().filter(|x| nodes.contains(x)).count()).unwrap_or(0));
    }
    let hubs: Vec<u128> = nodes.iter().copied().filter(|n| indeg[n] != 1 || outdeg[n] != 1).collect();
    let mut out = Vec::<SuperEdge>::new();
    for &u in &hubs {
        if let Some(vs0) = succ.get(&u) {
            for &v0 in vs0 {
                if !nodes.contains(&v0) { continue; }
                let mut chain = vec![u, v0];
                let mut cur = v0;
                while nodes.contains(&cur)
                    && indeg[&cur] == 1 && outdeg[&cur] == 1
                {
                    let next = succ.get(&cur).and_then(|v| v.first()).copied();
                    if let Some(nxt) = next {
                        if !nodes.contains(&nxt) { break; }
                        chain.push(nxt);
                        cur = nxt;
                        if indeg[&cur] != 1 || outdeg[&cur] != 1 { break; }
                    } else { break; }
                }
                out.push(SuperEdge { from: u, to: cur, chain });
            }
        }
    }
    (out, hubs)
}
fn index_from(edges: &[SuperEdge]) -> HashMap<u128, Vec<usize>> {
    let mut m: HashMap<u128, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        m.entry(e.from).or_default().push(i);
    }
    m
}
fn append_next_kmer(assembled: &mut String, next_seq: &str, k: usize) -> bool {
    if assembled.as_bytes()[(assembled.len() - (k - 1))..] == next_seq.as_bytes()[..(k - 1)] {
        assembled.push(next_seq.as_bytes()[k - 1] as char);
        true
    } else {
        let rc = revcomp(next_seq);
        if assembled.as_bytes()[(assembled.len() - (k - 1))..] == rc.as_bytes()[..(k - 1)] {
            assembled.push(rc.as_bytes()[k - 1] as char);
            true
        } else {
            false
        }
    }
}
fn assemble_chain_sequence(chain: &[u128], segs: &HashMap<u128, Segment>, k: usize) -> Option<String> {
    if chain.is_empty() { return None; }
    let mut seq = segs[&chain[0]].seq.clone();
    for node in chain.iter().skip(1) {
        if !append_next_kmer(&mut seq, &segs[node].seq, k) {
            // retry with RC of first
            let mut seq2 = revcomp(&segs[&chain[0]].seq);
            for node2 in chain.iter().skip(1) {
                if !append_next_kmer(&mut seq2, &segs[node2].seq, k) { return None; }
            }
            return Some(seq2);
        }
    }
    Some(seq)
}
fn enumerate_paths_and_sequences(
    start: u128,
    end: u128,
    edges: &[SuperEdge],
    segs: &HashMap<u128, Segment>,
    k: usize,
    max_paths: usize,
) -> Vec<String> {
    let from_map = index_from(edges);
    let mut out_seqs: Vec<String> = Vec::new();
    let mut stack: Vec<(u128, Vec<usize>, String)> = Vec::new();

    if let Some(eidxs) = from_map.get(&start) {
        for &ei in eidxs {
            if let Some(seq) = assemble_chain_sequence(&edges[ei].chain, segs, k) {
                stack.push((edges[ei].to, vec![ei], seq));
            }
        }
    }
    let guard = 2000usize;
    while let Some((at, idxs, seq)) = stack.pop() {
        if out_seqs.len() >= max_paths { break; }
        if at == end { out_seqs.push(seq.clone()); continue; }
        if idxs.len() > guard { continue; }
        if let Some(eidxs) = from_map.get(&at) {
            for &ei in eidxs {
                let e = &edges[ei];
                let mut seq2 = seq.clone();
                if assemble_chain_sequence(&e.chain, segs, k).is_some()
                    && append_next_kmer(&mut seq2, &segs[&e.chain[1]].seq, k)
                {
                    let mut idxs2 = idxs.clone();
                    idxs2.push(ei);
                    stack.push((e.to, idxs2, seq2));
                } else {
                    // slower but reliable: append whole chain
                    let mut ok = true;
                    for node in e.chain.iter().skip(1) {
                        if !append_next_kmer(&mut seq2, &segs[node].seq, k) { ok = false; break; }
                    }
                    if ok {
                        let mut idxs2 = idxs.clone();
                        idxs2.push(ei);
                        stack.push((e.to, idxs2, seq2));
                    }
                }
            }
        }
    }
    out_seqs
}

// BubbleGun shape
#[derive(Deserialize)]
struct BubbleChain { bubbles: Vec<Bubble> }
#[derive(Deserialize)]
struct Bubble { id: usize, ends: Vec<String>, inside: Vec<String> }

fn collect_bubbles_any(v: &Value) -> Vec<Bubble> {
    let mut out = Vec::new();
    match v {
        Value::Object(map) => {
            if let Some(b) = map.get("bubbles") {
                if let Ok(bc) = serde_json::from_value::<BubbleChain>(Value::Object(map.clone())) {
                    out.extend(bc.bubbles);
                } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(b.clone()) {
                    out.extend(bs);
                }
            } else {
                for vv in map.values() { out.extend(collect_bubbles_any(vv)); }
            }
        }
        Value::Array(arr) => {
            if let Ok(bcs) = serde_json::from_value::<Vec<BubbleChain>>(Value::Array(arr.clone())) {
                for bc in bcs { out.extend(bc.bubbles); }
            } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(Value::Array(arr.clone())) {
                out.extend(bs);
            } else {
                for vv in arr { out.extend(collect_bubbles_any(vv)); }
            }
        }
        _ => {}
    }
    out
}

// ====== High-k graph index (k-mer→nodes + adjacency) ======

fn base2(x: u8) -> Option<u8> {
    match x {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}
fn pack_kmer_u128(s: &str) -> Option<u128> {
    let mut v: u128 = 0;
    for &b in s.as_bytes() {
        let d = base2(b)?;
        v = (v << 2) | (d as u128);
    }
    Some(v)
}
fn roll_fwd(prev: u128, out_b: u8, k: usize) -> Option<u128> {
    let d = base2(out_b)?;
    let mask = if k == 64 { u128::MAX } else { (1u128 << (2*k)) - 1 };
    let v = ((prev << 2) | (d as u128)) & mask;
    Some(v)
}

struct HighIndex {
    k: usize,
    kmer2nodes: HashMap<u128, Vec<u128>>,
    succ: HashMap<u128, Vec<u128>>,
}
fn build_high_index(gfa: &str) -> HighIndex {
    let f = File::open(gfa).expect("open high-k GFA");
    let rdr = BufReader::new(f);
    let mut k: Option<usize> = None;
    let mut kmer2nodes: HashMap<u128, Vec<u128>> = HashMap::new();
    let mut succ: HashMap<u128, Vec<u128>> = HashMap::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                let id: u128 = cols[1].parse().unwrap();
                let seq = cols[2].as_bytes();
                if k.is_none() { k = Some(seq.len()); }
                if seq.len() != k.unwrap() {
                    // in a DBG S lines are kmers; if not, skip indexing this node
                    continue;
                }
                if let Some(code) = pack_kmer_u128(std::str::from_utf8(seq).unwrap()) {
                    kmer2nodes.entry(code).or_default().push(id);
                }
                succ.entry(id).or_default();
            }
            Some("L") => {
                let from: u128 = cols[1].parse().unwrap();
                let to: u128 = cols[3].parse().unwrap();
                succ.entry(from).or_default().push(to);
            }
            _ => {}
        }
    }
    HighIndex { k: k.expect("empty high-k graph"), kmer2nodes, succ }
}

// Walk candidate sequence through high-k graph (both strands); count solutions
fn count_walks(seq: &str, hi: &HighIndex, max_frontier: usize) -> usize {
    if seq.len() < hi.k { return 0; }
    let mut total = 0usize;
    for orient in 0..2 {
        let s = if orient == 0 { seq } else { &revcomp(seq) };
        let first = &s[..hi.k];
        let mut start_code = match pack_kmer_u128(first) { Some(c) => c, None => continue };
        if let Some(starts) = hi.kmer2nodes.get(&start_code) {
            // BFS frontier over nodes following sequence
            for &start_id in starts {
                let mut frontier: HashSet<u128> = HashSet::new();
                frontier.insert(start_id);
                let mut ok = true;
                for &b in s.as_bytes()[hi.k..].iter() {
                    let mut next_frontier: HashSet<u128> = HashSet::new();
                    // compute next kmer code from any representative (roll based on sequence itself)
                    start_code = roll_fwd(start_code, b, hi.k).unwrap();
                    if let Some(cands) = hi.kmer2nodes.get(&start_code) {
                        let cand_set: HashSet<u128> = cands.iter().copied().collect();
                        for &u in frontier.iter() {
                            if let Some(vs) = hi.succ.get(&u) {
                                for &v in vs {
                                    if cand_set.contains(&v) {
                                        next_frontier.insert(v);
                                        if next_frontier.len() >= max_frontier { break; }
                                    }
                                }
                            }
                            if next_frontier.len() >= max_frontier { break; }
                        }
                    }
                    if next_frontier.is_empty() { ok = false; break; }
                    frontier = next_frontier;
                }
                if ok { total += frontier.len().max(1); }
            }
        }
    }
    total
}

// ====== Main glue ======

fn main() {
    // Usage:
    // cargo run --release --bin project_highk_mapwalk -- \
    //   --low-gfa ecoli_21_graph.gfa \
    //   --bubbles ecoli_unlabeled_2.json \
    //   --high-gfa ecoli_41_graph.gfa \
    //   --out ecoli_dataset_missing_mapwalk.jsonl \
    //   [--allow-multi-occurrence] [--max-frontier 512] [--max-paths 1024]
    let mut low_gfa = None;
    let mut bubbles_json = None;
    let mut high_gfa = None;
    let mut out_path = None;
    let mut allow_multi = false;
    let mut max_frontier: usize = 512;
    let mut max_paths: usize = 1024;

    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--low-gfa" => low_gfa = args.next(),
            "--bubbles" => bubbles_json = args.next(),
            "--high-gfa" => high_gfa = args.next(),
            "--out" => out_path = args.next(),
            "--allow-multi-occurrence" => allow_multi = true,
            "--require-unique-occurrence" => allow_multi = false,
            "--max-frontier" => { max_frontier = args.next().unwrap().parse().unwrap(); }
            "--max-paths" => { max_paths = args.next().unwrap().parse().unwrap(); }
            x => { eprintln!("Unknown arg: {x}"); std::process::exit(2); }
        }
    }
    let low_gfa = low_gfa.expect("--low-gfa required");
    let bubbles_json = bubbles_json.expect("--bubbles required");
    let high_gfa = high_gfa.expect("--high-gfa required");
    let out_path = out_path.expect("--out required");

    eprintln!("[info] reading low-k graph: {low_gfa}");
    let (low_segs, low_links) = parse_gfa(&low_gfa);
    let (low_succ, low_pred) = build_adj(&low_segs, &low_links);
    let k_low = low_segs.values().next().expect("empty low-k").seq.len();

    eprintln!("[info] loading unlabeled BubbleGun file: {bubbles_json}");
    let bubbles_val: Value = serde_json::from_reader(File::open(&bubbles_json).expect("open bubbles")).expect("parse bubbles JSON");
    let bubbles = collect_bubbles_any(&bubbles_val);

    eprintln!("[info] building high-k index from: {high_gfa}");
    let hi = build_high_index(&high_gfa);
    eprintln!("[info] high-k k={}", hi.k);

    let mut out = std::fs::File::create(&out_path).expect("create out");
    let mut labeled = 0usize;

    for b in bubbles {
        if b.ends.len() != 2 { continue; }
        let start = parse_u128(&b.ends[0]);
        let end   = parse_u128(&b.ends[1]);

        // induced nodes: inside + endpoints
        let mut nodes: HashSet<u128> = b.inside.iter().map(|s| parse_u128(s)).collect();
        nodes.insert(start);
        nodes.insert(end);

        let (superedges, _supernodes) = compact_with_chains(&nodes, &low_succ, &low_pred);
        let seqs = enumerate_paths_and_sequences(start, end, &superedges, &low_segs, k_low, max_paths);

        if seqs.is_empty() { continue; }

        // count walks per candidate
        let mut counts: Vec<usize> = Vec::with_capacity(seqs.len());
        for s in seqs.iter() {
            counts.push(count_walks(s, &hi, max_frontier));
        }
        // winner logic
        let mut winner: Option<(usize, usize)> = None; // (idx, hits)
        let mut clash = false;
        for (i, &h) in counts.iter().enumerate() {
            if h > 0 {
                if winner.is_none() { winner = Some((i, h)); }
                else { clash = true; break; }
            }
        }
        if let Some((i, hits)) = winner {
            if !clash && (allow_multi || hits == 1) {
                let reason = if hits == 1 { "highk_mapwalk_unique" } else { "highk_mapwalk_multi_occ" };
                let obj = serde_json::json!({
                    "bubble_id": b.id,
                    "label_path": i,
                    "label_reason": reason,
                    "winner_occurrences": hits
                });
                writeln!(out, "{}", obj.to_string()).unwrap();
                labeled += 1;
            }
        }
    }

    eprintln!("[done] labeled {labeled} bubbles -> {out_path}");
}
