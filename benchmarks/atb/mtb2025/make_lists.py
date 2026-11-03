#!/usr/bin/env python3
"""Sample lists for the Mtb 2024-08 / 2025-05 build from ATB assemblies.

Groups: K = Mtb in both releases (2024-08 species_calls 'Species', 2025-05 'sylph_species'),
D = Mtb in 2024-08 only (dropped), N = Mtb in 2025-05 only (new).
Writes lists/samples.tsv (sample, group, archive, path in archive), lists/archives.tsv (archive, url,
md5, size_MB, n_K, n_D, n_N, n_other) and lists/<group>.txt.
"""
import collections, csv, gzip, os
# The ATB metadata: 2024-08/ and 2025-05/species_calls.tsv.gz, files/file_list.all.latest.tsv.gz
M = os.environ["ATB_META"]
W = os.path.dirname(os.path.abspath(__file__))
MTB = "Mycobacterium tuberculosis"
new = {r["sample_accession"] for r in csv.DictReader(gzip.open(f"{M}/2025-05/species_calls.tsv.gz", "rt"), delimiter="\t") if r["sylph_species"] == MTB}
old = {r["Sample"] for r in csv.DictReader(gzip.open(f"{M}/2024-08/species_calls.tsv.gz", "rt"), delimiter="\t") if r["Species"] == MTB}
group = {s: ("K" if s in new else "D") for s in old}
group.update({s: "N" for s in new - old})
arch = {}; count = collections.defaultdict(collections.Counter); rows = []
for r in csv.DictReader(gzip.open(f"{M}/files/file_list.all.latest.tsv.gz", "rt"), delimiter="\t"):
    a = r["tar_xz"]
    arch[a] = (r["tar_xz_url"], r["tar_xz_md5"], r["tar_xz_size_MB"])
    g = group.get(r["sample"])
    count[a][g or "other"] += 1
    if g:
        rows.append((r["sample"], g, a, r["filename_in_tar_xz"]))
assert len(rows) == len(group), (len(rows), len(group))
os.makedirs(f"{W}/lists", exist_ok=True)
with open(f"{W}/lists/samples.tsv", "w") as o:
    o.write("sample\tgroup\tarchive\tpath\n")
    for r in sorted(rows):
        o.write("\t".join(r) + "\n")
with open(f"{W}/lists/archives.tsv", "w") as o:
    o.write("archive\turl\tmd5\tsize_MB\tn_K\tn_D\tn_N\tn_other\n")
    for a in sorted(a for a in count if sum(count[a][g] for g in "KDN")):
        c = count[a]
        o.write(f"{a}\t{arch[a][0]}\t{arch[a][1]}\t{arch[a][2]}\t{c['K']}\t{c['D']}\t{c['N']}\t{c['other']}\n")
for g in "KDN":
    with open(f"{W}/lists/{g}.txt", "w") as o:
        o.writelines(s + "\n" for s, gg, _, _ in sorted(rows) if gg == g)
print({g: sum(1 for r in rows if r[1] == g) for g in "KDN"}, "archives:", sum(1 for a in count if sum(count[a][g] for g in "KDN")))
