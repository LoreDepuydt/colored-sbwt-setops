#!/bin/bash
# The benchmark of intersect and difference. As queue_merge.bash: every run through run.bash,
# skipped if it is in the ledger already, repetitions round by round (E. coli REPS_ECOLI, pneumo REPS_PNEUMO).
# The inputs are the indexes made by queue_merge.bash (indexes/<dataset>/merge/) and the full indexes.
#
# For each pair (A, B), three runs: intersect --colors union --merge-shared-colors (<pair>_union),
# intersect --colors intersect (<pair>_inter), difference A - B (<pair>_diff).
#   E. coli pairs: full with 1 genome (one), 10 % (tenth), 50 % (half), all but 1 genome (allbut1);
#     two indexes that share 10 % or 50 % of the genomes (shared10, shared50); two disjoint halves
#     (disjoint, no --colors intersect run: no shared colors); half and full (reverse: half - full is empty)
#   S. pneumoniae pairs: full with 1 genome (one), full with 50 % (half)
# There is no direct build to compare with: the results are not the index of a set of genomes.
# Checks, in round 1 (checks.tsv):
#   - B a part of A (the pairs with the full index), E. coli: intersect --colors intersect = B exactly
#   - every pair: the k-mers add up, |A - B| + |A ∩ B| = |A|, and both color modes keep the same k-mers
#     (themisto2 stats of the results, not timed; the counts go to setops_stats.tsv)
# Caps: E. coli 16 GiB, pneumo 240 GiB. Usage: bash queue_setops.bash   (settings: THREADS (32),
#   REPS_ECOLI (5), REPS_PNEUMO (3); for a small test: PNEUMO=0, ECOLI_FULL, EI, ECOLI_PAIRS, W_ROOT)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS_ECOLI=${REPS_ECOLI:-5}
REPS_PNEUMO=${REPS_PNEUMO:-3}
PNEUMO=${PNEUMO:-1}
CHECKS=${CHECKS:-$P/checks.tsv}
STATS=${STATS:-$P/setops_stats.tsv}
source "$P/lib.bash"
W_ROOT=${W_ROOT:-$B/work}
ECOLI_FULL=${ECOLI_FULL:-$B/indexes/ecoli/scale_3682.thm2}
EI=${EI:-$B/indexes/ecoli/merge}
PNEUMO_FULL=$B/indexes/pneumo/s_pneumoniae.thm2
PI=$B/indexes/pneumo/merge
# <name>:<A>:<B>, with A and B as paths
ECOLI_PAIRS=${ECOLI_PAIRS:-"one:$ECOLI_FULL:$EI/disjoint_1_B.thm2 tenth:$ECOLI_FULL:$EI/disjoint_368_B.thm2 half:$ECOLI_FULL:$EI/parts2_0.thm2 allbut1:$ECOLI_FULL:$EI/disjoint_1_A.thm2 shared10:$EI/shared_10_A.thm2:$EI/shared_10_B.thm2 shared50:$EI/shared_50_A.thm2:$EI/shared_50_B.thm2 disjoint:$EI/parts2_0.thm2:$EI/parts2_1.thm2 reverse:$EI/parts2_0.thm2:$ECOLI_FULL"}
PNEUMO_PAIRS="one:$PNEUMO_FULL:$PI/rest_121033.thm2 half:$PNEUMO_FULL:$PI/first_60517.thm2"
[[ -s $STATS ]] || printf 'index\tn_kmers\tn_colors\tn_distinct_color_sets\tfull_index_size_bytes\tdate\n' > "$STATS"

# kmers <index> <cap>: number of k-mers of an index (themisto2 stats, capped, not timed; kept in STATS)
kmers() {
    local n
    n=$(awk -F'\t' -v i="$1" '$1 == i { n = $2 } END { print n }' "$STATS")
    if [[ -z $n ]]; then
        systemd-run --user --scope --quiet -p MemoryMax="$2G" -p MemorySwapMax=0 nice -n 19 \
            "$BIN/themisto2" stats -i "$1" -t "$THREADS" 2> /dev/null > "$W/stats.out" || return 1
        awk -F'\t' -v i="$1" -v d="$(date -Is)" '{ v[$1] = $2 } END { printf "%s\t%s\t%s\t%s\t%s\t%s\n", i, v["n_kmers"], v["n_colors"], v["n_distinct_color_sets"], v["full_index_size_bytes"], d }' "$W/stats.out" >> "$STATS"
        n=$(awk -F'\t' '$1 == "n_kmers" { print $2 }' "$W/stats.out")
    fi
    echo "$n"
}
# record <description> <result>: one line in checks.tsv (a count check: no color names column)
record() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "n/a" "$(date -Is)" | tee -a "$CHECKS"
}
done_check() { [[ -f $CHECKS ]] && awk -F'\t' -v d="$1" '$1 == d { f = 1 } END { exit !f }' "$CHECKS"; }

# pair <dataset> <name> <A> <B> <round> <cap>: the three runs of a pair; in round 1 the results are
# kept for the checks, then deleted
pair() {
    local ds=$1 nm=$2 a=$3 b=$4 r=$5 cap=$6 u=$W/$2.union.thm2 i=$W/$2.inter.thm2 d=$W/$2.diff.thm2
    CAP_GB=$cap INPUTS="$a $b" run "${ds}_intersect_${nm}_union_r$r" "$u" 0 \
        "$BIN/themisto2" intersect --index1 "$a" --index2 "$b" --colors union --merge-shared-colors -o "$u" -t "$THREADS" --no-lcs || log "FAILED ${ds}_intersect_${nm}_union_r$r"
    # Disjoint pair: no shared colors, so --colors intersect has nothing to keep (and its result,
    # with 0 colors, does not load); counted with the union result instead
    if [[ $nm == disjoint ]]; then
        i=$u
    else
        CAP_GB=$cap INPUTS="$a $b" run "${ds}_intersect_${nm}_inter_r$r" "$i" 0 \
            "$BIN/themisto2" intersect --index1 "$a" --index2 "$b" --colors intersect -o "$i" -t "$THREADS" --no-lcs || log "FAILED ${ds}_intersect_${nm}_inter_r$r"
    fi
    CAP_GB=$cap INPUTS="$a $b" run "${ds}_difference_${nm}_r$r" "$d" 0 \
        "$BIN/themisto2" difference --index1 "$a" --index2 "$b" -o "$d" -t "$THREADS" --no-lcs || log "FAILED ${ds}_difference_${nm}_r$r"
    if (( r == 1 )); then
        local desc="$ds $nm: |A - B| + |A ∩ B| = |A|, union and intersect modes the same k-mers"
        if ! done_check "$desc"; then
            if [[ -f $u && -f $i && -f $d ]]; then
                local na nu ni nd
                na=$(kmers "$a" "$cap") && nu=$(kmers "$u" "$cap") && ni=$(kmers "$i" "$cap") && nd=$(kmers "$d" "$cap") &&
                    record "$desc" "$( (( nd + ni == na && nu == ni )) && echo equal || echo "DIFFERENT ($nd + $ni vs $na; union $nu)")" ||
                    log "no count check of $ds $nm: stats failed"
            else
                log "no count check of $ds $nm: a result is missing"
            fi
        fi
        # B is a part of A (the pairs with the full index first): intersect --colors intersect = B
        [[ $ds == ecoli && $a == "$ECOLI_FULL" ]] && CAP_GB=16 compare_once "$ds intersect $nm --colors intersect = index2" "$i" "$b"
    fi
    rm -f -- "$u" "$i" "$d"
}
compare_once() {
    done_check "$1" && return 0
    [[ -f $2 && -f $3 ]] || { log "no check of '$1': $2 or $3 missing"; return 1; }
    check "$CHECKS" "$1" "$2" "$3"
}

round() {
    local ds=$1 r=$2 cap=$3 p nm a b
    W=$W_ROOT/${ds}_setops; mkdir -p "$W/tmp"
    for p in $4; do
        IFS=: read -r nm a b <<< "$p"
        [[ -f $a && -f $b ]] || { log "FAILED $ds $nm: input $a or $b missing (run queue_merge.bash first)"; continue; }
        pair "$ds" "$nm" "$a" "$b" "$r" "$cap"
    done
}

for (( r = 1; r <= REPS_ECOLI || (PNEUMO && r <= REPS_PNEUMO); r++ )); do
    log "round $r"
    (( r <= REPS_ECOLI )) && round ecoli "$r" 16 "$ECOLI_PAIRS"
    (( PNEUMO && r <= REPS_PNEUMO )) && round pneumo "$r" 240 "$PNEUMO_PAIRS"
done
log "queue_setops done"
