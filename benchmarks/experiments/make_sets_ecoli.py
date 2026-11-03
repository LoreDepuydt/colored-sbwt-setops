#!/usr/bin/env python3
"""The genome sets of the E. coli benchmark, as lists of color names (genome accessions), one per line,
into sets/ecoli/. All sets are cut from one random order of the 3,682 genomes (seed 1), so they are
reproducible and nested where that matters (a prefix of a random order is a random set).

  order.txt                 all genomes, in the random order
  scale_<N>                 the first N: the indexes of the scale series (500, 1000, 2000, 3682)
  keep_<n>                  subset keep sets of the full index: the first n (n genomes, or n % of them)
  merge_disjoint_<m>_A/_B   A = all but the last m, B = the last m (B is the update)
  merge_shared_<pct>_A/_B   two halves that share <pct> % of all genomes (A = first half + overlap,
                            B = second half + overlap), A union B = all
  parts<K>_<i>              all genomes in K equal consecutive parts (multi-way and incremental merge)
  nested_<n>                the first n, as index2 of difference/intersect against the full index

Usage: python3 make_sets_ecoli.py <folder of the genomes (*.fna)> <out folder>
"""
import os
import random
import sys

genome_dir, out = sys.argv[1], sys.argv[2]
names = sorted(f[:-len(".fna")] for f in os.listdir(genome_dir) if f.endswith(".fna"))
assert len(names) == 3682, len(names)
order = names[:]
random.Random(1).shuffle(order)
N = len(order)
os.makedirs(out, exist_ok=True)

def write(name, items):
    assert len(set(items)) == len(items), name
    with open(os.path.join(out, name + ".txt"), "w") as f:
        f.writelines(x + "\n" for x in items)

write("order", order)
for n in (500, 1000, 2000, N):
    write(f"scale_{n}", order[:n])
for n in (1, 10, 100, N // 10, N // 2, N - N // 10):
    write(f"keep_{n}", order[:n])
for m in (N // 2, N // 10, N // 100, 1):
    write(f"merge_disjoint_{m}_A", order[:N - m])
    write(f"merge_disjoint_{m}_B", order[N - m:])
h = N // 2
for pct in (10, 50, 100):
    s = round(pct / 100 * N / 2)          # each half gets s genomes of the other one
    a, b = order[:min(N, h + s)], order[max(0, h - s):]
    write(f"merge_shared_{pct}_A", a)
    write(f"merge_shared_{pct}_B", b)
    assert set(a) | set(b) == set(order)
for k in (2, 4, 8, 16):
    bounds = [round(i * N / k) for i in range(k + 1)]
    for i in range(k):
        write(f"parts{k}_{i}", order[bounds[i]:bounds[i + 1]])
for n in (1, N // 10, N // 2):
    write(f"nested_{n}", order[:n])
print(f"{len(os.listdir(out))} sets in {out}")
