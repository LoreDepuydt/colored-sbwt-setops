#!/usr/bin/env python3
"""Helpers of mtb_case.bash (small, single-threaded; standard library only).

  percolour <prefix> <names.tsv> <out.tsv>        canonical k-mers per colour of an export (GFA + colour sets),
                                                  and per number of colours in the set
  species   <percolour.tsv> <out.tsv>             per ATB species: shared k-mers with Mtb, taxonomy, ATB samples,
                                                  Mtb samples in which sylph co-detects it; sums per category;
                                                  Spearman of shared k-mers against co-detections
  simreads  <ref.fa> <n> <len> <seed> <out.fa>    error-free reads, both strands
  fpsummary <label> <pseudoalignment> <n_reads>   reads with >= 1 colour

GFA segments are canonical unitigs, so k-mer counts here are canonical (half of an index's n_kmers).
"""
import collections, csv, gzip, os, random, re, sys

K = 31
# The ATB metadata: atb_species_counts.2024-08.tsv, gtdb/bac120_taxonomy_r214.tsv.gz,
# 2024-08/species_calls.tsv.gz and 2024-08/sylph.tsv.gz
META = os.environ["ATB_META"]
MTB = "Mycobacterium tuberculosis"


def read_names(path):
    names = {}
    for line in open(path):
        i, n = line.rstrip("\n").split("\t", 1)
        names[int(i)] = n
    return names


def colour_species(name):
    m = re.match(r"per_species_unitigs/(.+)-unitigs-k31\.fna$", name)
    return m.group(1) if m else name


def percolour(prefix, names_path, out):
    names = read_names(names_path)
    sets = {}
    for line in open(prefix + ".color_sets.txt"):
        f = line.split()
        sets[int(f[0].split("=")[1])] = [int(x) for x in f[2:]]
    per_set = collections.Counter()
    for line in open(prefix + ".unitigs.gfa"):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            per_set[next(int(x[5:]) for x in f[3:] if x.startswith("cs:i:"))] += len(f[2]) - K + 1
    per_colour, per_size = collections.Counter(), collections.Counter()
    for s, km in per_set.items():
        per_size[min(len(sets[s]), 11)] += km
        for c in sets[s]:
            per_colour[c] += km
    total = sum(per_set.values())
    with open(out, "w") as o:
        o.write(f"# canonical k-mers: {total}\nkind\tname\tcanonical_kmers\tshare\n")
        for n, km in sorted(per_size.items()):
            o.write(f"set_size\t{n if n < 11 else '>10'}\t{km}\t{km / total:.6f}\n")
        for c, km in per_colour.most_common():
            o.write(f"colour\t{colour_species(names[c])}\t{km}\t{km / total:.6f}\n")
    print(f"{out}: {total} canonical k-mers in {len(per_colour)} colours; colours per k-mer: "
          + ", ".join(f"{n if n < 11 else '>10'}: {100 * km / total:.1f}%" for n, km in sorted(per_size.items())))


def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for t in range(i, j + 1):
            r[order[t]] = (i + j) / 2
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx, vy = sum((a - mx) ** 2 for a in rx) ** 0.5, sum((b - my) ** 2 for b in ry) ** 0.5
    return cov / (vx * vy) if vx and vy else float("nan")


def species(percolour_path, out):
    shared = {}
    total = None
    for line in open(percolour_path):
        if line.startswith("# canonical k-mers:"):
            total = int(line.split(":")[1])
        f = line.rstrip("\n").split("\t")
        if f[0] == "colour":
            shared[f[1]] = int(f[2])
    col = {r["colour"]: r for r in csv.DictReader(open(f"{META}/atb_species_counts.2024-08.tsv"), delimiter="\t")}
    tax = {}
    with gzip.open(f"{META}/gtdb/bac120_taxonomy_r214.tsv.gz", "rt") as f:
        for line in f:
            t = dict(x.split("__", 1) for x in line.rstrip("\n").split("\t")[1].split(";"))
            tax[t["s"]] = t
    mtb_samples = set()
    with gzip.open(f"{META}/2024-08/species_calls.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Species"] == MTB:
                mtb_samples.add(r["Sample"])
    co = collections.defaultdict(set)
    with gzip.open(f"{META}/2024-08/sylph.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Sample"] in mtb_samples and r["Species"] != MTB:
                co[r["Species"]].add(r["Sample"])
    rows = []
    for c, r in col.items():
        sp = r["species"]
        if sp == MTB:
            continue
        t = tax.get(sp, {})
        n_co = len(co.get(sp, ()))
        cat = ("Mycobacteriaceae" if t.get("f") == "Mycobacteriaceae" else
               "mixed_call_with_Mtb" if ";" in sp and MTB in sp.split(";") else
               "codetected_ge50" if n_co >= 50 else "codetected_1_49" if n_co > 0 else "never_codetected")
        rows.append((sp, c, shared.get(c, 0), int(r["n_samples"]), n_co, cat, t.get("o", ""), t.get("f", ""), t.get("g", "")))
    rows.sort(key=lambda r: -r[2])
    with open(out, "w") as o:
        o.write("species\tcolour\tshared_canonical_kmers_with_mtb\tn_atb_samples\tn_mtb_samples_codetected\tcategory\torder\tfamily\tgenus\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    print(f"{out}: {sum(1 for r in rows if r[2])} species share k-mers with the Mtb colour (O has {total} canonical k-mers)")
    by_cat = collections.Counter()
    for r in rows:
        by_cat[r[5]] += r[2]
    print("  shared k-mers summed per category (a k-mer counts once per species):",
          ", ".join(f"{c} {k}" for c, k in by_cat.most_common()))
    print("  top 15 species:")
    for r in rows[:15]:
        print(f"    {r[2]:>12d}  co-detected in {r[4]:>5d} Mtb samples  {r[5]:<20} {r[0]}")
    x = [r[2] for r in rows]
    y = [r[4] for r in rows]
    print(f"  Spearman, shared k-mers vs sylph co-detections, all {len(rows)} species: {spearman(x, y):.3f}")
    sub = [r for r in rows if r[4] > 0]
    if len(sub) > 2:
        print(f"  ... over the {len(sub)} co-detected species: {spearman([r[2] for r in sub], [r[4] for r in sub]):.3f}")


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def simreads(ref, n, length, seed, out):
    seqs = []
    for block in open(ref).read().split(">")[1:]:
        s = re.sub(r"[^ACGT]", "N", "".join(block.split("\n")[1:]).upper())
        if len(s) >= int(length):
            seqs.append(s)
    tot = sum(len(s) for s in seqs)
    rnd = random.Random(int(seed))
    with open(out, "w") as o:
        i = 0
        while i < int(n):
            x = rnd.randrange(tot)
            for s in seqs:
                if x < len(s):
                    break
                x -= len(s)
            r = s[x:x + int(length)]
            if len(r) < int(length) or "N" in r:
                continue
            if rnd.random() < 0.5:
                r = revcomp(r)
            o.write(f">r{i}\n{r}\n"); i += 1


def fpsummary(label, path, n_reads):
    hit = sum(1 for line in open(path) if len(line.split()) > 1)
    print(f"{label}\t{hit}\t{n_reads}\t{100 * hit / int(n_reads):.3f}")


if __name__ == "__main__":
    {"percolour": percolour, "species": species, "simreads": simreads, "fpsummary": fpsummary}[sys.argv[1]](*sys.argv[2:])
