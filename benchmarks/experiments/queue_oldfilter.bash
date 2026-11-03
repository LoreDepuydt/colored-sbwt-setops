#!/bin/bash
# The repair of the SBWT of our earlier work (an auxiliary SBWT and a three-way interleaving: the binary
# of tag recomb27-threeway, OLD) against the current one (insertion), on the same inputs and with the
# same commands as the current runs. Every run through run.bash, skipped if in the ledger, round by round; the outputs are
# deleted after each run (their md5 is in the ledger, to compare with the current runs).
#   (a) subsets of queue_subset.bash: E. coli compl_1 (remove 1),
#       compl_100 (remove 100), keep_1841 (50 %), keep_368 (10 %), keep_1; S. pneumoniae compl_1, keep_12103
#       (10 %). Runs <ds>_subsetold_<set>_r<i>.
#   (b) the E. coli pairs of queue_setops.bash, the same three runs per pair (no --colors intersect for
#       disjoint). Runs ecoli_intersectold_<pair>_{union,inter}_r<i>, ecoli_differenceold_<pair>_r<i>.
#   (c) S. pneumoniae difference full - 1 genome (queue_setops.bash "one", difference only), cap 240 GiB
#       as the current run. Runs pneumo_differenceold_one_r<i>.
# Order: (c) first, all its rounds; then the E. coli rounds; the S. pneumoniae subsets of (a) last.
# Rounds: E. coli REPS_ECOLI (5), S. pneumoniae REPS_PNEUMO (3). Caps as the current runs: E. coli 16 GiB,
# pneumo 160 GiB. Summary: python3 oldfilter_stats.py.
# Usage: bash queue_oldfilter.bash   (settings: OLD ($BIN/oldfilter/themisto2), THREADS (32), REPS_ECOLI,
#   REPS_PNEUMO, PNEUMO (1); DRY=1 prints the runs only)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS_ECOLI=${REPS_ECOLI:-5}
REPS_PNEUMO=${REPS_PNEUMO:-3}
PNEUMO=${PNEUMO:-1}
DRY=${DRY:-0}
source "$P/lib.bash"
OLD=${OLD:-$BIN/oldfilter/themisto2}
[[ -x $OLD ]] || { log "ABORT: no binary $OLD (build tag recomb27-threeway)"; exit 1; }
log "earlier repair: $OLD, md5 $(md5sum < "$OLD" | cut -d' ' -f1)"
W=$B/work/oldfilter_runs; mkdir -p "$W/tmp"
ECOLI_FULL=$B/indexes/ecoli/scale_3682.thm2
EI=$B/indexes/ecoli/merge
PNEUMO_FULL=$B/indexes/pneumo/s_pneumoniae.thm2
ECOLI_SUBSETS="compl_1 compl_100 keep_1841 keep_368 keep_1"
PNEUMO_SUBSETS="compl_1 keep_12103"
PI=$B/indexes/pneumo/merge
ECOLI_PAIRS="one:$ECOLI_FULL:$EI/disjoint_1_B.thm2 tenth:$ECOLI_FULL:$EI/disjoint_368_B.thm2 half:$ECOLI_FULL:$EI/parts2_0.thm2 allbut1:$ECOLI_FULL:$EI/disjoint_1_A.thm2 shared10:$EI/shared_10_A.thm2:$EI/shared_10_B.thm2 shared50:$EI/shared_50_A.thm2:$EI/shared_50_B.thm2 disjoint:$EI/parts2_0.thm2:$EI/parts2_1.thm2 reverse:$EI/parts2_0.thm2:$ECOLI_FULL"

# go <name> <cap> <inputs> <output> <args...>: one timed run of the old binary, output deleted
go() {
    local name=$1 cap=$2 in=$3 out=$4; shift 4
    if (( DRY )); then done_run "$name" && echo "skip $name" || echo "$name (cap $cap): themisto2-old $*"; return 0; fi
    CAP_GB=$cap INPUTS="$in" run "$name" "$out" 1 "$OLD" "$@" -o "$out" -t "$THREADS" --no-lcs || log "FAILED $name"
}
# subset_flags <dataset> <set>: the flags of queue_subset.bash (compl_<n> = --complement of keep_<n>)
subset_flags() {
    local n=${2#*_}
    [[ $2 == compl_* ]] && echo "--complement --colors-file $P/sets/$1/keep_$n.txt" || echo "--colors-file $P/sets/$1/keep_$n.txt"
}

for f in "$ECOLI_FULL" "$PNEUMO_FULL" "$PI/rest_121033.thm2"; do [[ -f $f ]] || { log "ABORT: no $f"; exit 1; }; done
for p in $ECOLI_PAIRS; do IFS=: read -r _ a b <<< "$p"; [[ -f $a && -f $b ]] || { log "ABORT: missing input of $p"; exit 1; }; done
if (( PNEUMO )); then
    for (( r = 1; r <= REPS_PNEUMO; r++ )); do
        go "pneumo_differenceold_one_r$r" 240 "$PNEUMO_FULL $PI/rest_121033.thm2" "$W/d.thm2" difference --index1 "$PNEUMO_FULL" --index2 "$PI/rest_121033.thm2"
    done
fi
for (( r = 1; r <= REPS_ECOLI; r++ )); do
    log "E. coli round $r"
    for s in $ECOLI_SUBSETS; do
        # shellcheck disable=SC2046
        go "ecoli_subsetold_${s}_r$r" 16 "$ECOLI_FULL" "$W/sub.thm2" subset --index "$ECOLI_FULL" $(subset_flags ecoli "$s")
    done
    for p in $ECOLI_PAIRS; do
        IFS=: read -r nm a b <<< "$p"
        go "ecoli_intersectold_${nm}_union_r$r" 16 "$a $b" "$W/u.thm2" intersect --index1 "$a" --index2 "$b" --colors union --merge-shared-colors
        [[ $nm == disjoint ]] || go "ecoli_intersectold_${nm}_inter_r$r" 16 "$a $b" "$W/i.thm2" intersect --index1 "$a" --index2 "$b" --colors intersect
        go "ecoli_differenceold_${nm}_r$r" 16 "$a $b" "$W/d.thm2" difference --index1 "$a" --index2 "$b"
    done
done
if (( PNEUMO )); then
    for (( r = 1; r <= REPS_PNEUMO; r++ )); do
        for s in $PNEUMO_SUBSETS; do
            # shellcheck disable=SC2046
            go "pneumo_subsetold_${s}_r$r" 160 "$PNEUMO_FULL" "$W/sub.thm2" subset --index "$PNEUMO_FULL" $(subset_flags pneumo "$s")
        done
    done
fi
log "queue_oldfilter done"
