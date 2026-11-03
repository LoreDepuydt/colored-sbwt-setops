#!/usr/bin/env python3
"""The genome sets of the S. pneumoniae benchmark, as lists of color names of the index (e.g.
s_pneumoniae/SAMEA1024559.fa), one per line, into sets/pneumo/. As make_sets_ecoli.py: one random order
of the 121,034 colors (seed 1), every set cut from it (nested).

  order.txt       all colors, in the random order
  keep_<n>        the first n: subset keep sets, and the genomes of the direct builds (1, 10, 100,
                  1,000, 10 % = 12,103, 50 %, 90 %)
  atb.tsv         per color: color name, ATB accession, size and md5 (S3 ETag) of <accession>.fa.gz in
                  the ATB assemblies bucket (from the inventory of 2025-02-16)

Usage: python3 make_sets_pneumo.py <color_names.tsv (id TAB name)> <ATB inventory .csv.gz> <out folder>
"""
import csv
import gzip
import os
import random
import sys

names_file, inventory, out = sys.argv[1:4]
with open(names_file) as f:
    names = [line.rstrip("\n").split("\t")[1] for line in f if line.strip()]
assert len(names) == 121034 and len(set(names)) == len(names), len(names)
acc = {n: os.path.basename(n).removesuffix(".fa") for n in names}
wanted = set(acc.values())
atb = {}
with gzip.open(inventory, "rt") as f:
    for row in csv.reader(f):
        a = row[1].removesuffix(".fa.gz")
        if a in wanted:
            atb[a] = (row[2], row[4])
missing = wanted - atb.keys()
assert not missing, f"{len(missing)} colors not in ATB, e.g. {sorted(missing)[:3]}"

order = names[:]
random.Random(1).shuffle(order)
N = len(order)
os.makedirs(out, exist_ok=True)
def write(name, items):
    with open(os.path.join(out, name + ".txt"), "w") as f:
        f.writelines(x + "\n" for x in items)
write("order", order)
for n in (1, 10, 100, 1000, round(N / 10), round(N / 2), round(N * 9 / 10)):
    write(f"keep_{n}", order[:n])
with open(os.path.join(out, "atb.tsv"), "w") as f:
    f.write("color\taccession\tbytes\tmd5\n")
    for n in order:
        f.write(f"{n}\t{acc[n]}\t{atb[acc[n]][0]}\t{atb[acc[n]][1]}\n")
print(f"{len(os.listdir(out))} files in {out}")
