#!/usr/bin/env python3
"""The earlier repair of the SBWT (queue_oldfilter.bash, runs *old_*) against the current one (the runs of
queue_subset.bash and queue_setops.bash with the same inputs), from ledger.tsv and the logs in runs/.

One row per operation (subsets and E. coli intersect/difference), medians over the repetitions:
  wall_s, peak_GiB      the ledger (time -v wall time, max RSS)
  sbwt_s                the SBWT step of the operation: "Pass 1: computing incoming-edge coverage" to
                        "Building subset rank structure" (the filter of a subset, or the intersection/
                        difference of the SBWTs), with the dummy repair
  repair_s              the dummy repair alone: "<n> source(s) missing dummy chains" to "Dummy repair:
                        building ... bit-rows via three-way pass" (0 when no source misses its chain)
  sources, orphans, |Q| the sources missing dummy chains, the orphans of the dead-chain removal that get new
                        chains (the same in both filters), and their sum
  inserted              current filter: "inserting <n> new dummy node(s)"; the old filter has no such line
  aux_nodes             old filter: the nodes of its auxiliary SBWT (from the interleaving with it)
  same_md5              every repetition of both filters wrote the same index (md5 of the ledger)
  n1, n2                the positions the operation works on, as the log gives them: for intersect and
                        difference the nodes (SBWT sets) of the two inputs ("Computing merge interleaving:
                        index1 has <n1> sets ..., index2 has <n2> sets"), whose merge plan has at most
                        n1 + n2 positions; for a subset the nodes of the filter's input, the full index (its
                        sbwt_length; the subset log has no such line), and n2 empty
  result_sets           the nodes of the result SBWT ("Number of sets in <op> SBWT: <n>")
  mode                  intersect: the color mode (--colors intersect, the default, or --colors union
                        --merge-shared-colors); empty for difference and subset
  cur_runs              the dates of the current runs the numbers come from (no current run was redone)
Writes oldfilter_compare.tsv. Usage: python3 oldfilter_stats.py   (in paper/)
"""
import csv, re, statistics, sys
import os, socket
from collections import defaultdict
from pathlib import Path

P = Path(__file__).resolve().parent
# The host of the runs, as run.bash names their folders (runs/<name>_<host>)
HOST = os.environ.get("BENCH_HOST") or socket.gethostname().split(".")[0]
OLD = re.compile(r"(ecoli|pneumo)_(subset|intersect|difference)old_(.+)_r(\d+)$")


# The nodes of the full indexes (sbwt_length), the input of the subsets
FULL_SETS = {"ecoli": 342587824, "pneumo": 417856502}


def parse(log):
    lines = []
    for line in log.read_text(errors="replace").splitlines():
        t, _, msg = line.partition(" ")
        try:
            lines.append((float(t), msg))
        except ValueError:
            pass
    first = lambda pat: next(((t, m) for t, m in lines if re.search(pat, m)), (None, None))
    v = {}
    a, _ = first(r"Pass 1: computing incoming-edge coverage")
    b, _ = first(r"Building subset rank structure")
    v["sbwt_s"] = b - a if a is not None and b is not None else None
    t0, m = first(r"(\d+) source\(s\) missing dummy chains")
    t1, _ = first(r"Dummy repair: building \S+ SBWT bit-rows via three-way pass")
    v["sources"] = int(re.search(r"(\d+) source\(s\) missing", m).group(1)) if m else 0
    v["repair_s"] = t1 - t0 if t0 is not None and t1 is not None else 0.0
    _, m = first(r"(\d+) orphan\(s\)")
    v["orphans"] = int(re.search(r"(\d+) orphan\(s\)", m).group(1)) if m else 0
    _, m = first(r"inserting (\d+) new dummy node\(s\)")
    v["inserted"] = int(re.search(r"inserting (\d+)", m).group(1)) if m else (0 if not v["sources"] and not v["orphans"] else None)
    # the auxiliary SBWT of the old filter: the second index of the interleaving after it is built
    aux = None
    after = [m for t, m in lines if t0 is not None and t >= t0]
    for m in after:
        x = re.search(r"aux=(\d+) sets", m) or re.search(r"Computing merge interleaving: index1 has \d+ sets \(\d+ k-mers\), index2 has (\d+) sets", m)
        if x:
            aux = int(x.group(1))
            break
    v["aux_nodes"] = aux
    _, m = first(r"Computing merge interleaving: index1 has \d+ sets")
    x = m and re.search(r"index1 has (\d+) sets.*index2 has (\d+) sets", m)
    v["n1"], v["n2"] = (int(x.group(1)), int(x.group(2))) if x else (None, None)
    _, m = first(r"Number of sets in \S+ SBWT: \d+")
    v["result_sets"] = int(re.search(r"SBWT: (\d+)", m).group(1)) if m else None
    return v


def main():
    led = {}
    with open(P / "ledger.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["host"] == HOST and r["exit"] == "0" and r["killed"] in ("", "0"):
                led[r["name"]] = r
    groups = defaultdict(lambda: {"old": [], "cur": []})
    for name in led:
        m = OLD.match(name)
        if m:
            ds, op, rest, rnd = m.groups()
            key = f"{ds}_{op}_{rest}"
            groups[key]["old"].append(name)
            cur = f"{key}_r{rnd}"
    for key in groups:
        groups[key]["cur"] = sorted(n for n in led if re.fullmatch(re.escape(key) + r"_r\d+", n))
    rows, problems = [], []
    med = lambda xs: round(statistics.median(xs), 2) if xs and None not in xs else ""
    for key in sorted(groups):
        row = {"operation": key}
        md5s = set()
        for side in ("cur", "old"):
            names = groups[key][side]
            parsed = []
            for n in names:
                try:
                    parsed.append(parse(P / "runs" / f"{n}_{HOST}" / "stderr.log"))
                except OSError as e:
                    problems.append(f"{n}: {e}")
            md5s |= {led[n]["output_md5"] for n in names}
            row[f"n_{side}"] = len(names)
            row[f"wall_s_{side}"] = med([float(led[n]["wall_s"]) for n in names])
            row[f"peak_GiB_{side}"] = med([int(led[n]["max_rss_kb"]) / 2**20 for n in names])
            row[f"sbwt_s_{side}"] = med([p["sbwt_s"] for p in parsed])
            row[f"repair_s_{side}"] = med([p["repair_s"] for p in parsed])
            if "_subset_" not in key:
                for c in ("n1", "n2", "result_sets"):
                    vals = {p[c] for p in parsed}
                    if len(vals) == 1 and None not in vals:
                        row[c] = vals.pop()
                    elif vals:
                        problems.append(f"{key} {side}: {c} {vals}")
            else:
                vals = {p["result_sets"] for p in parsed}
                if len(vals) == 1 and None not in vals:
                    row["result_sets"] = vals.pop()
            for c in ("sources", "orphans"):
                vals = {p[c] for p in parsed}
                if len(vals) > 1:
                    problems.append(f"{key} {side}: {c} differ between repetitions: {vals}")
                row[f"{c}_{side}"] = vals.pop() if len(vals) == 1 else ""
            if side == "cur":
                vals = {p["inserted"] for p in parsed}
                row["inserted_cur"] = vals.pop() if len(vals) == 1 else ""
            else:
                vals = {p["aux_nodes"] for p in parsed}
                row["aux_nodes_old"] = vals.pop() if len(vals) == 1 else ""
        for c in ("sources", "orphans"):
            if row[f"{c}_cur"] != row[f"{c}_old"] and row[f"n_old"] and row[f"n_cur"]:
                problems.append(f"{key}: {c} {row[f'{c}_cur']} (current) vs {row[f'{c}_old']} (old)")
        q = [row.get("sources_cur"), row.get("orphans_cur")]
        row["Q"] = sum(q) if all(isinstance(x, int) for x in q) else ""
        row["same_md5"] = "yes" if len(md5s) == 1 and row["n_old"] and row["n_cur"] else ("NO" if len(md5s) > 1 else "")
        if row["same_md5"] == "NO":
            problems.append(f"{key}: different indexes ({len(md5s)} md5s)")
        ds = key.split("_")[0]
        if "_subset_" in key:
            row["n1"], row["n2"] = FULL_SETS[ds], ""
        row["mode"] = "--colors intersect" if key.endswith("_inter") else ("--colors union --merge-shared-colors" if key.endswith("_union") else "")
        row["cur_runs"] = ", ".join(sorted({led[n]["start"][:10] for n in groups[key]["cur"]}))
        rows.append(row)
    if not rows:
        print("no old-filter runs in the ledger yet")
        return 0
    cols = ["operation", "mode", "n1", "n2", "result_sets", "cur_runs", "n_cur", "n_old", "wall_s_cur", "wall_s_old", "peak_GiB_cur", "peak_GiB_old", "sbwt_s_cur", "sbwt_s_old",
            "repair_s_cur", "repair_s_old", "sources_cur", "orphans_cur", "Q", "inserted_cur", "aux_nodes_old", "same_md5"]
    with open(P / "oldfilter_compare.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} operations; wrote oldfilter_compare.tsv")
    for p in problems:
        print("PROBLEM:", p, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
