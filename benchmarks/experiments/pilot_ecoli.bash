#!/bin/bash
# Pilot of the E. coli benchmark: tests run.bash, the pipelines and the correctness checks on one run
# of each step, before the full queue. Every run goes through run.bash (wait for a quiet server,
# cap, samples, ledger; cold start: its inputs are dropped from the page cache first), named pilot_*. Runs already in the ledger are skipped, so it can be restarted.
#   1. the full index (3,682 genomes) with the direct pipeline: ggcat -> sbwt build -> themisto2 build
#   2. the two halves parts2_0 and parts2_1 by subset of the full index
#   3. merge of the halves             check: = the full index
#   4. intersect full, parts2_0 with --colors intersect      check: = parts2_0
#   5. difference full - parts2_0      (no exact reference here; the full queue compares to a direct build)
#   6. subset keep_368, with and without --low-ram-mode       check: the two are the same
# Checks: export of both indexes, compare_color_dumps (the k-mers and the color set of every k-mer),
# and the color names; one line per check in pilot_checks.tsv.
# Usage: bash pilot_ecoli.bash   (settings: THREADS (32), CAP_GB (16))
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
GENOMES=${BENCH_DIR:-/tmp/colored-sbwt-bench}/data/ecoli/genomes
IDX=${BENCH_DIR:-/tmp/colored-sbwt-bench}/indexes/ecoli          # input indexes: kept, never deleted by the runs
W=${BENCH_DIR:-/tmp/colored-sbwt-bench}/work/pilot               # results and temporary files of the runs: emptied after each run
SETS=$P/sets/ecoli
THREADS=${THREADS:-32}
CAP_GB=${CAP_GB:-16}
CHECKS=$P/pilot_checks.tsv
source "$P/lib.bash"
mkdir -p "$IDX" "$W/tmp" "$B/work/fof"
# fof <set>: genome<TAB>name lines of a set, in the order of the set
fof() { awk -v g="$GENOMES" '{ print g "/" $1 ".fna\t" $1 }' "$SETS/$1.txt" > "$B/work/fof/$1.tsv"; echo "$B/work/fof/$1.tsv"; }

# 1. Full index, README pipeline
full=$IDX/scale_3682.thm2
F=$(fof scale_3682)
cut -f1 "$F" > "$B/work/fof/scale_3682.genomes.txt"
INPUTS="@$F" run pilot_ggcat_3682 "$W/full.unitigs.fa" 0 "$BIN/ggcat" build -l "$B/work/fof/scale_3682.genomes.txt" -s 1 -k $K -m 8 -p -j "$THREADS" -t "$W/tmp" -o "$W/full.unitigs.fa" || exit 1
INPUTS="$W/full.unitigs.fa" run pilot_sbwt_3682 "$W/full.sbwt" 0 "$BIN/sbwt" build -i "$W/full.unitigs.fa" -o "$W/full" -k $K -t "$THREADS" -m 8 -r -v --temp-dir "$W/tmp" || exit 1
INPUTS="@$F $W/full.sbwt" run pilot_build_3682 "$full" 0 "$BIN/themisto2" build --file-colors "$F" -s "$W/full.sbwt" -o "$full" -k $K -d 30 -t "$THREADS" --no-lcs || exit 1
rm -f -- "$W/full.unitigs.fa" "$W/full.sbwt"

# 2. Halves by subset
for h in parts2_0 parts2_1; do
    INPUTS="$full" run "pilot_subset_$h" "$IDX/$h.thm2" 0 "$BIN/themisto2" subset --index "$full" --colors-file "$SETS/$h.txt" -o "$IDX/$h.thm2" -t "$THREADS" --no-lcs || exit 1
done

# 3. Merge
printf '%s\n' "$IDX/parts2_0.thm2" "$IDX/parts2_1.thm2" > "$B/work/fof/merge_parts2.txt"
INPUTS="$IDX/parts2_0.thm2 $IDX/parts2_1.thm2" run pilot_merge_parts2 "$W/merge.thm2" 0 "$BIN/themisto2" merge --index-file-list "$B/work/fof/merge_parts2.txt" --temp-dir "$W/tmp" -o "$W/merge.thm2" -t "$THREADS" --no-lcs &&
    check "$CHECKS" "merge(parts2_0, parts2_1) = full" "$W/merge.thm2" "$full"
rm -f -- "$W/merge.thm2"

# 4. Intersect with --colors intersect
INPUTS="$full $IDX/parts2_0.thm2" run pilot_intersect_full_parts2_0 "$W/inter.thm2" 0 "$BIN/themisto2" intersect --index1 "$full" --index2 "$IDX/parts2_0.thm2" --colors intersect -o "$W/inter.thm2" -t "$THREADS" --no-lcs &&
    check "$CHECKS" "intersect(full, parts2_0, --colors intersect) = parts2_0" "$W/inter.thm2" "$IDX/parts2_0.thm2"
rm -f -- "$W/inter.thm2"

# 5. Difference
INPUTS="$full $IDX/parts2_0.thm2" run pilot_difference_full_parts2_0 "$W/diff.thm2" 1 "$BIN/themisto2" difference --index1 "$full" --index2 "$IDX/parts2_0.thm2" -o "$W/diff.thm2" -t "$THREADS" --no-lcs

# 6. Subset, default and --low-ram-mode
INPUTS="$full" run pilot_subset_keep_368 "$W/s368.thm2" 0 "$BIN/themisto2" subset --index "$full" --colors-file "$SETS/keep_368.txt" -o "$W/s368.thm2" -t "$THREADS" --no-lcs &&
    INPUTS="$full" run pilot_subset_keep_368_lowram "$W/s368_lowram.thm2" 0 "$BIN/themisto2" subset --index "$full" --colors-file "$SETS/keep_368.txt" -o "$W/s368_lowram.thm2" -t "$THREADS" --no-lcs --low-ram-mode &&
    check "$CHECKS" "subset keep_368 = subset keep_368 --low-ram-mode" "$W/s368.thm2" "$W/s368_lowram.thm2"
rm -f -- "$W/s368.thm2" "$W/s368_lowram.thm2"
log "pilot done"
