// Unitig length x GC histogram, tips and simple bubbles of a themisto2 `export --gfa` file.
// The graph is bidirected: one segment per canonical unitig, links are k-1 overlaps between
// segment sides (left = start of the forward sequence, right = its end).
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

const NONE: u32 = u32::MAX;

fn parse_id(s: &[u8]) -> u64 {
    // "u123"
    let mut x = 0u64;
    for &c in &s[1..] {
        x = x * 10 + (c - b'0') as u64;
    }
    x
}

fn len_class(kmers: u32, k: u32) -> &'static str {
    match kmers {
        1 => "a:1",
        2..=10 => "b:2-10",
        11..=30 => "c:11-30",
        _ if kmers == k => "d:31(=k)",
        32..=61 => "e:32-61",
        62..=100 => "f:62-100",
        101..=300 => "g:101-300",
        301..=1000 => "h:301-1000",
        1001..=10000 => "i:1001-10000",
        _ => "j:>10000",
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let gfa = &args[1];
    let out = &args[2];
    let k: u32 = 31;
    let overlap = (k - 1) as usize;

    let mut seg_id: Vec<u64> = Vec::new();
    let mut seg_kmers: Vec<u32> = Vec::new();
    let mut seg_gc: Vec<u16> = Vec::new(); // GC per mille
    let mut links: Vec<(u64, bool, u64, bool)> = Vec::new(); // (a, a is +, b, b is +)
    let mut non_plus = 0u64;

    let mut r = BufReader::with_capacity(1 << 24, File::open(gfa).unwrap());
    let mut line = Vec::new();
    let mut n_lines = 0u64;
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        n_lines += 1;
        if n_lines % 50_000_000 == 0 {
            eprintln!("{} lines, {} segments, {} links", n_lines, seg_id.len(), links.len());
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        let mut f = line.split(|&c| c == b'\t');
        match f.next() {
            Some(b"S") => {
                let id = parse_id(f.next().unwrap());
                let seq = f.next().unwrap();
                let gc = seq.iter().filter(|&&c| c == b'G' || c == b'C' || c == b'g' || c == b'c').count();
                seg_id.push(id);
                seg_kmers.push((seq.len() - overlap) as u32);
                seg_gc.push(((gc as u64 * 1000) / seq.len() as u64) as u16);
            }
            Some(b"L") => {
                let a = parse_id(f.next().unwrap());
                let oa = f.next().unwrap();
                let b = parse_id(f.next().unwrap());
                let ob = f.next().unwrap();
                if oa != b"+" || ob != b"+" {
                    non_plus += 1;
                }
                links.push((a, oa == b"+", b, ob == b"+"));
            }
            _ => {}
        }
    }
    let n = seg_id.len();
    eprintln!("read {} segments, {} links", n, links.len());

    // Sort segments by id to map ids to indexes.
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by_key(|&i| seg_id[i as usize]);
    let sorted_ids: Vec<u64> = order.iter().map(|&i| seg_id[i as usize]).collect();
    drop(seg_id);
    let idx = |id: u64| -> u32 { order[sorted_ids.binary_search(&id).expect("link to unknown segment")] };

    // Sides: 2i = left end of segment i, 2i+1 = right end. "L a oa b ob" leaves a by its right
    // end if oa is +, else by its left end, and enters b by its left end if ob is +, else by its
    // right end. Each adjacency may be listed in both directions, so the pairs are deduplicated.
    let mut edges: Vec<(u32, u32)> = links
        .iter()
        .map(|&(a, pa, b, pb)| {
            let sa = 2 * idx(a) + if pa { 1 } else { 0 };
            let sb = 2 * idx(b) + if pb { 0 } else { 1 };
            (sa.min(sb), sa.max(sb))
        })
        .collect();
    let n_links = links.len();
    drop(links);
    edges.sort_unstable();
    edges.dedup();
    let n_adjacencies = edges.len();
    let mut deg = vec![0u32; 2 * n];
    let mut nbr = vec![NONE; 2 * n]; // the neighbour side, meaningful when deg == 1
    let mut self_loops = 0u64;
    for &(s, t) in &edges {
        if s / 2 == t / 2 {
            self_loops += 1;
        }
        deg[s as usize] += 1;
        nbr[s as usize] = t;
        if s != t {
            deg[t as usize] += 1;
            nbr[t as usize] = s;
        }
    }
    drop(edges);

    // Categories.
    // 0 isolated (both sides free), 1 tip (exactly one side free), 2 bubble branch, 3 other.
    let mut cat = vec![3u8; n];
    for i in 0..n {
        let (dl, dr) = (deg[2 * i], deg[2 * i + 1]);
        if dl == 0 && dr == 0 {
            cat[i] = 0;
        } else if dl == 0 || dr == 0 {
            cat[i] = 1;
        }
    }

    // Simple bubbles: >= 2 segments with one neighbour on each side and the same pair of
    // neighbour sides (unordered, since branches can be stored in either orientation).
    let mut cand: Vec<(u32, u32, u32)> = Vec::new(); // (neighbour side, neighbour side, segment)
    for i in 0..n {
        if deg[2 * i] == 1 && deg[2 * i + 1] == 1 {
            let (x, y) = (nbr[2 * i], nbr[2 * i + 1]);
            if x != y && x / 2 != i as u32 && y / 2 != i as u32 {
                cand.push((x.min(y), x.max(y), i as u32));
            }
        }
    }
    cand.sort_unstable();
    let mut n_bubbles = 0u64;
    let mut n_bubbles_snp = 0u64; // two branches, both exactly k k-mers
    let mut n_bubbles_two = 0u64;
    let mut n_bubbles_indel = 0u64; // two branches of different lengths
    let mut branch_count_hist: HashMap<usize, u64> = HashMap::new();
    let mut j = 0;
    while j < cand.len() {
        let mut e = j + 1;
        while e < cand.len() && cand[e].0 == cand[j].0 && cand[e].1 == cand[j].1 {
            e += 1;
        }
        let m = e - j;
        if m >= 2 {
            n_bubbles += 1;
            *branch_count_hist.entry(m).or_insert(0) += 1;
            for t in j..e {
                cat[cand[t].2 as usize] = 2;
            }
            if m == 2 {
                n_bubbles_two += 1;
                let (l1, l2) = (seg_kmers[cand[j].2 as usize], seg_kmers[cand[j + 1].2 as usize]);
                if l1 == k && l2 == k {
                    n_bubbles_snp += 1;
                } else if l1 != l2 {
                    n_bubbles_indel += 1;
                }
            }
        }
        j = e;
    }
    drop(cand);

    // Histograms.
    let cat_name = ["isolated", "tip", "bubble_branch", "other"];
    let mut lg: HashMap<(&str, u16), (u64, u64)> = HashMap::new(); // (len class, GC 5% bin)
    let mut cg: HashMap<(u8, u16), (u64, u64)> = HashMap::new(); // (category, GC 1% bin)
    let mut cl: HashMap<(u8, &str), (u64, u64)> = HashMap::new(); // (category, len class)
    let mut tot = [(0u64, 0u64); 4];
    for i in 0..n {
        let km = seg_kmers[i] as u64;
        let lc = len_class(seg_kmers[i], k);
        let e = lg.entry((lc, seg_gc[i] / 50)).or_insert((0, 0));
        e.0 += 1;
        e.1 += km;
        let e = cg.entry((cat[i], seg_gc[i] / 10)).or_insert((0, 0));
        e.0 += 1;
        e.1 += km;
        let e = cl.entry((cat[i], lc)).or_insert((0, 0));
        e.0 += 1;
        e.1 += km;
        tot[cat[i] as usize].0 += 1;
        tot[cat[i] as usize].1 += km;
    }

    let mut w = BufWriter::new(File::create(format!("{}.length_gc.tsv", out)).unwrap());
    writeln!(w, "length_class_kmers\tgc_bin_pct\tn_unitigs\tn_kmers").unwrap();
    let mut v: Vec<_> = lg.into_iter().collect();
    v.sort();
    for ((lc, b), (u, km)) in v {
        writeln!(w, "{}\t{}-{}\t{}\t{}", lc, b * 5, b * 5 + 5, u, km).unwrap();
    }
    let mut w = BufWriter::new(File::create(format!("{}.category_gc.tsv", out)).unwrap());
    writeln!(w, "category\tgc_pct\tn_unitigs\tn_kmers").unwrap();
    let mut v: Vec<_> = cg.into_iter().collect();
    v.sort();
    for ((c, b), (u, km)) in v {
        writeln!(w, "{}\t{}\t{}\t{}", cat_name[c as usize], b, u, km).unwrap();
    }
    let mut w = BufWriter::new(File::create(format!("{}.category_length.tsv", out)).unwrap());
    writeln!(w, "category\tlength_class_kmers\tn_unitigs\tn_kmers").unwrap();
    let mut v: Vec<_> = cl.into_iter().collect();
    v.sort();
    for ((c, lc), (u, km)) in v {
        writeln!(w, "{}\t{}\t{}\t{}", cat_name[c as usize], lc, u, km).unwrap();
    }

    let mut w = BufWriter::new(File::create(format!("{}.summary.txt", out)).unwrap());
    let total_kmers: u64 = seg_kmers.iter().map(|&x| x as u64).sum();
    writeln!(w, "segments\t{}\nlinks\t{}\nnon_plus_links\t{}\nadjacencies\t{}\nself_loops\t{}\ntotal_kmers\t{}", n, n_links, non_plus, n_adjacencies, self_loops, total_kmers).unwrap();
    for c in 0..4 {
        writeln!(w, "{}\tunitigs {}\tkmers {}", cat_name[c], tot[c].0, tot[c].1).unwrap();
    }
    writeln!(w, "bubbles\t{}\ntwo_branch\t{}\ntwo_branch_both_k_kmers(snp-like)\t{}\ntwo_branch_unequal(indel-like)\t{}", n_bubbles, n_bubbles_two, n_bubbles_snp, n_bubbles_indel).unwrap();
    let mut bh: Vec<_> = branch_count_hist.into_iter().collect();
    bh.sort();
    for (m, c) in bh {
        writeln!(w, "bubbles_with_{}_branches\t{}", m, c).unwrap();
    }
    let mut dh: HashMap<(u32, u32), u64> = HashMap::new();
    for i in 0..n {
        *dh.entry((deg[2 * i].min(deg[2 * i + 1]).min(5), deg[2 * i].max(deg[2 * i + 1]).min(5))).or_insert(0) += 1;
    }
    let mut dv: Vec<_> = dh.into_iter().collect();
    dv.sort();
    for ((a, b), c) in dv {
        writeln!(w, "side_degrees_{}_{}\t{}", a, b, c).unwrap();
    }
}
