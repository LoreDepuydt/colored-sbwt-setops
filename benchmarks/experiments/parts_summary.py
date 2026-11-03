#!/usr/bin/env python3
"""Build times of the merged parts (queue_parts.bash) next to the random sets of the same size
(queue_subset.bash). From ledger.tsv (host HOST, exit 0). For each repetition of a direct build, the
time is the sum of its three steps (ggcat, sbwt, themisto) and the peak the largest step peak; the table
gives the medians over the repetitions, and the faster of direct and plain build (plain only for at most
100 genomes). Writes parts_build.tsv.
Usage: python3 parts_summary.py   (in paper/)
"""
import csv, re, statistics
import os, socket
from collections import defaultdict
from pathlib import Path

P = Path(__file__).resolve().parent
# The host of the runs, as run.bash names their folders (runs/<name>_<host>)
HOST = os.environ.get("BENCH_HOST") or socket.gethostname().split(".")[0]
ROWS = [  # (row, part runs, the random set used so far)
    ("E. coli add 1", "ecoli_part{}_merge_disjoint_1_B", "ecoli_{}_keep_1"),
    ("E. coli add 368", "ecoli_part{}_merge_disjoint_368_B", "ecoli_{}_keep_368"),
    ("S. pneumoniae add 1", "pneumo_part{}_rest_121033", "pneumo_{}_keep_1"),
]

runs = {}
with open(P / "ledger.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["host"] == HOST and r["exit"] == "0" and r["killed"] in ("", "0"):
            runs[r["name"]] = (float(r["wall_s"]), int(r["max_rss_kb"]) / 2**20)


def direct(pre):
    """median sum of the steps and median of the largest step peak, over the rounds with all 3 steps"""
    by_round = defaultdict(dict)
    for name, v in runs.items():
        m = re.fullmatch(re.escape(pre) + r"_(ggcat|sbwt|themisto)_r(\d+)", name)
        if m:
            by_round[m.group(2)][m.group(1)] = v
    full = [d for d in by_round.values() if len(d) == 3]
    if not full:
        return None
    return (len(full), statistics.median(sum(w for w, _ in d.values()) for d in full),
            statistics.median(max(p for _, p in d.values()) for d in full))


def plain(pre):
    v = [runs[n] for n in runs if re.fullmatch(re.escape(pre) + r"_r\d+", n)]
    if not v:
        return None
    return len(v), statistics.median(w for w, _ in v), statistics.median(p for _, p in v)


def fmt(x):
    return "" if x is None else f"{x[1]:.2f} s, {x[2]:.2f} GiB (n={x[0]})"


def best(d, p):
    c = [(x, k) for x, k in ((d, "direct"), (p, "plain")) if x is not None]
    return min(c, key=lambda t: t[0][1]) if c else (None, "")


out = []
for row, part, rand in ROWS:
    pd, pp = direct(part.format("direct")), plain(part.format("plain"))
    rd, rp = direct(rand.format("direct")), plain(rand.format("plainbuild"))
    (bp, kp), (br, kr) = best(pd, pp), best(rd, rp)
    out.append({"row": row, "part_direct": fmt(pd), "part_plain": fmt(pp), "part_best": f"{fmt(bp)} {kp}".strip(),
                "random_direct": fmt(rd), "random_plain": fmt(rp), "random_best": f"{fmt(br)} {kr}".strip()})
with open(P / "parts_build.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
for r in out:
    print(r["row"], "| part:", r["part_best"] or "(not run yet)", "| random set:", r["random_best"])
