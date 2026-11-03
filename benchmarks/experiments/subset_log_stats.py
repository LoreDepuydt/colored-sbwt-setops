#!/usr/bin/env python3
"""Dummy repair counts and time per phase of the subset runs, read from their logs (no new runs).

For every ecoli_subset_* / pneumo_subset_* run in ledger.tsv (host HOST, exit 0, not killed), reads
runs/<name>_<HOST>/stderr.log, where every line starts with the epoch at which it was written.
Writes, one row per run without the _r<i> suffix:
  repair_stats.tsv  the counts of sbwt::set_operations::intersect_difference and "Keeping <n> of <N>
                    k-mers" (the same in every repetition, else the run is reported and left out),
                    and the key k-mers of PHASE 1/3 ("Marked <x> % of all k-mers")
  phase_times.tsv   the median over the repetitions of the time of each phase (s):
    load        start of the command to "Index loaded" (cold start: the index is read from disk)
    dbg_marks   to "Marking the k-mers with the kept colors"
    mark_walk   to "Keeping ... k-mers" (marking the k-mers with the kept colors, walking the unitigs)
    filter      to "=== Phase 1/3" (filtering the SBWT, with the dummy repair, and the DBG of the subset)
      repair    the part of filter from "... source(s) missing dummy chains" to the three-way pass
                (which replaces Pass 2 of the filter without repair); 0 when no source misses its chain
    phase1      to "=== PHASE 2/3" (marking the key k-mers)
    phase2      to "=== PHASE 3/3" (color set fingerprints)
    phase3      to "Serializing subset" (distinct color set storage)
    serialize   to the end (start + wall_s of the ledger; 0 if that is up to 50 ms early, see parse)
The start is the epoch of the first line (the cgroup line, written when the command starts).
Usage: python3 subset_log_stats.py   (in paper/)
"""
import csv, re, statistics, sys
import os, socket
from collections import defaultdict
from pathlib import Path

P = Path(__file__).resolve().parent
# The host of the runs, as run.bash names their folders (runs/<name>_<host>)
HOST = os.environ.get("BENCH_HOST") or socket.gethostname().split(".")[0]
NUM = r"(\d+)"
COUNTS = {  # column: (pattern, groups); a count is 0 when its line is absent (no repair)
    "sources_missing_chains": (r"\[filter\] " + NUM + r" source\(s\) missing dummy chains", 1),
    "dead_end_dummies": (r"Dead-end dummies \(no surviving outgoing edge\): " + NUM, 1),
    "dead_chain_nodes_marked": (r"Dead-chain removal: " + NUM + r" node\(s\) marked from", 1),
    "dead_chain_dead_ends": (r"marked from " + NUM + r" dead end\(s\)", 1),
    "orphans_new_chains": (r"dead end\(s\); " + NUM + r" orphan\(s\) get new dummy chains", 1),
    "inserted_dummies": (r"Dummy repair: inserting " + NUM + r" new dummy node\(s\)", 1),
    "kept_kmers": (r"Keeping " + NUM + r" of \d+ k-mers", 1),
    "total_kmers": (r"Keeping \d+ of " + NUM + r" k-mers", 1),
    "key_kmers_pct": (r"Marked ([\d.]+) % of all k-mers", 1),
}
MARKS = [  # phase: the line that ends it
    ("load", "Index loaded"),
    ("dbg_marks", "Marking the k-mers with the kept colors"),
    ("mark_walk", re.compile(r"Keeping \d+ of \d+ k-mers")),
    ("filter", "=== Phase 1/3"),
    ("phase1", "=== PHASE 2/3"),
    ("phase2", "=== PHASE 3/3"),
    ("phase3", "Serializing subset"),
]
REPAIR = ("source(s) missing dummy chains", "Dummy repair: building filter SBWT bit-rows via three-way pass")
PHASES = [m[0] for m in MARKS] + ["serialize"]


def parse(log, wall):
    lines = []
    for line in log.read_text().splitlines():
        t, _, msg = line.partition(" ")
        lines.append((float(t), msg))
    if not lines[0][1].startswith("cgroup "):
        raise ValueError("the first line is not the cgroup line")
    text = "\n".join(m for _, m in lines)
    counts = {}
    for col, (pat, g) in COUNTS.items():
        found = re.findall(pat, text)
        if len(found) > 1:
            raise ValueError(f"{col}: {len(found)} matches")
        counts[col] = found[0] if found else ("" if col == "key_kmers_pct" else "0")
    t0 = lines[0][0]
    ends, i = [], 0
    for _, mark in MARKS:
        hit = (lambda m: mark.search(m)) if isinstance(mark, re.Pattern) else (lambda m: mark in m)
        while i < len(lines) and not hit(lines[i][1]):
            i += 1
        if i == len(lines):
            raise ValueError(f"no '{mark if isinstance(mark, str) else mark.pattern}' after the previous phase")
        ends.append(lines[i][0])
    # wall_s is rounded to 10 ms, and a tiny result is written in a few ms: a few ms early is 0 s
    if t0 + wall < ends[-1] - 0.05:
        raise ValueError("start + wall_s is before the serializing line")
    ends.append(max(t0 + wall, ends[-1]))
    times = dict(zip(PHASES, [b - a for a, b in zip([t0] + ends, ends)]))
    rs = [t for t, m in lines if REPAIR[0] in m]
    re_ = [t for t, m in lines if REPAIR[1] in m]
    times["repair"] = (re_[0] - rs[0]) if rs and re_ else 0.0
    if rs and not re_:
        raise ValueError("repair start without three-way pass")
    return counts, times


def main():
    runs = {}
    with open(P / "ledger.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if re.match(r"(ecoli|pneumo)_subset_", r["name"]) and r["host"] == HOST and r["exit"] == "0" \
                    and r["killed"] in ("", "0"):
                runs[r["name"]] = float(r["wall_s"])  # the last line of a name wins
    groups = defaultdict(list)
    for name in runs:
        groups[re.sub(r"_r\d+$", "", name)].append(name)
    problems = []
    rep_rows, phase_rows = [], []
    for g in sorted(groups):
        parsed = []
        for name in sorted(groups[g]):
            try:
                parsed.append((name, *parse(P / "runs" / f"{name}_{HOST}" / "stderr.log", runs[name])))
            except (OSError, ValueError, IndexError) as e:
                problems.append(f"{name}: {e}")
        if not parsed:
            continue
        variants = {tuple(c.items()) for _, c, _ in parsed}
        if len(variants) > 1:
            problems.append(f"{g}: the counts differ between repetitions: " +
                            "; ".join(f"{n} {dict(c)}" for n, c, _ in parsed))
        rep_rows.append({"run": g, "n_reps": len(parsed), "same_in_all_reps": "yes" if len(variants) == 1 else "NO",
                         **parsed[0][1]})
        med = {p: statistics.median(t[p] for _, _, t in parsed) for p in PHASES + ["repair"]}
        phase_rows.append({"run": g, "n_reps": len(parsed), "wall_s": round(statistics.median(runs[n] for n, _, _ in parsed), 2),
                           **{p: round(med[p], 2) for p in PHASES[:4]}, "repair": round(med["repair"], 2),
                           **{p: round(med[p], 2) for p in PHASES[4:]}})
    for out, rows in (("repair_stats.tsv", rep_rows), ("phase_times.tsv", phase_rows)):
        with open(P / out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    print(f"{len(runs)} runs in {len(groups)} groups; wrote repair_stats.tsv and phase_times.tsv")
    for p in problems:
        print("PROBLEM:", p, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
