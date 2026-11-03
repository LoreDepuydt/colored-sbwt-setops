#!/bin/bash
# Are the ATB assemblies the genomes of s_pneumoniae.thm2? Builds an index directly from the downloaded
# genomes of a keep set (ATB .fa.gz files, with the color names of the index), takes the subset of
# s_pneumoniae.thm2 to the same colors, and compares the two: the same k-mers and color set of every
# k-mer means the downloads are the genomes the index was built from. One line in pneumo_checks.tsv.
# Both runs through run.bash (pneumo_check_*), skipped if in the ledger. The subset loads the 43 GB
# index: needs CAP_GB (60) + 24 GiB available and a quiet server.
# Usage: bash pneumo_check.bash [keep set, default keep_10]   (first: bash download_pneumo.bash <set>)
#   BUILD_ONLY=1: only the direct build (small), not the subset
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
SET=${1:-keep_10}
THREADS=${THREADS:-32}
CAP_GB=${CAP_GB:-60}
W=${BENCH_DIR:-/tmp/colored-sbwt-bench}/work/pneumo_check/$SET
source "$P/lib.bash"
INDEX=$B/indexes/pneumo/s_pneumoniae.thm2
GENOMES=$B/data/pneumo/genomes
mkdir -p "$W/tmp" "$B/work/fof"

# The genomes with the color names of the index, in the order of the index (as the subset numbers them)
F=$B/work/fof/pneumo_$SET.tsv
awk -F'\t' -v g="$GENOMES" 'NR == FNR { s[$1] = 1; next } ($2 in s) { a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' \
    "$P/sets/pneumo/$SET.txt" "$P/sets/pneumo/color_names.tsv" > "$F"
missing=$(cut -f1 "$F" | while read -r f; do [[ -f $f ]] || echo "$f"; done | wc -l)
(( missing == 0 && $(wc -l < "$F") == $(wc -l < "$P/sets/pneumo/$SET.txt") )) || { log "$missing genomes missing: bash download_pneumo.bash $SET"; exit 1; }

CAP_GB=16 run "pneumo_check_build_$SET" "$W/direct.thm2" 0 \
    "$BIN/themisto2" build --file-colors "$F" -o "$W/direct.thm2" -k $K -t "$THREADS" --temp-dir "$W/tmp" --no-lcs || exit 1
[[ ${BUILD_ONLY:-0} == 1 ]] && { log "BUILD_ONLY=1: stop before the subset"; exit 0; }
[[ -f $INDEX ]] || { log "no index at $INDEX"; exit 1; }
run "pneumo_check_subset_$SET" "$W/subset.thm2" 0 \
    "$BIN/themisto2" subset --index "$INDEX" --colors-file "$P/sets/pneumo/$SET.txt" -o "$W/subset.thm2" -t "$THREADS" --no-lcs || exit 1
CAP_GB=16 check "$P/pneumo_checks.tsv" "ATB genomes of $SET = subset of s_pneumoniae.thm2" "$W/direct.thm2" "$W/subset.thm2"
