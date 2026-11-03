#!/bin/bash
# The merge benchmark: merge, against building the union directly; plus plain themisto2 build for the
# small sets of queue_subset.bash. As queue_subset.bash: every run through run.bash, skipped
# if it is in the ledger already, repetitions round by round (E. coli REPS_ECOLI, pneumo REPS_PNEUMO).
#
# Inputs (made once, by subset of the full index, runs named *_input_*, kept in indexes/):
#   E. coli: the A and B indexes of sets/ecoli/merge_* and the parts of sets/ecoli/parts<K>_<i>
#   S. pneumoniae: the first n colors of s_pneumoniae.thm2 in its own order, and the others
# Every merge here has all genomes of its dataset as its union, in the order of the full index, so its
# result must equal the full index (check in lib.bash, in round 1, one line in checks.tsv). S. pneumoniae
# results are checked on a subset: both the result and the full index are cut to the 1,000 genomes of
# sets/pneumo/keep_1000.txt (not timed), then compared; a whole S. pneumoniae export is over 800 GB.
#
# E. coli, each round:
#   disjoint: merge A + B, B = the last m genomes, m = 1 (all + 1), 36 (99/1), 368 (90/10), 1841 (50/50)
#   shared colors: A and B share 10, 50 or 100 % of the genomes (--merge-shared-colors)
#   multi-way: the K parts in one call, K = 2, 4, 8, 16
#   one by one: the K parts added one at a time (K - 1 merges, each a run _s<i>), K = 4, 8, 16
#   50/50 also with --keep-redundant-dummies and with --low-ram-mode
#   the direct build of the union (= the full index; rounds 2-5, round 1 is the build of the input)
#   plain themisto2 build of keep_1, keep_10, keep_100 (the fastest build for so few genomes)
# S. pneumoniae, each round: disjoint merges of the first n colors and the others, n = 121033 (all + 1),
#   60517 (50/50); plain themisto2 build of keep_1, keep_10, keep_100. No direct build of the union
#   (its genomes are not downloaded).
# Caps: E. coli 16 GiB; pneumo subsets (inputs) 160 GiB, merges MERGE_CAP_PNEUMO (240) GiB, builds 16 GiB.
#
# Usage: bash queue_merge.bash   (settings: THREADS (32), REPS_ECOLI (5), REPS_PNEUMO (3); for a small
#   test: ECOLI_FULL, SETS, PNEUMO=0 (E. coli only), DISJOINT, SHARED, PARTS, W_ROOT, IDX_ROOT)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS_ECOLI=${REPS_ECOLI:-5}
REPS_PNEUMO=${REPS_PNEUMO:-3}
PNEUMO=${PNEUMO:-1}
MERGE_CAP_PNEUMO=${MERGE_CAP_PNEUMO:-240}
CHECKS=${CHECKS:-$P/checks.tsv}
source "$P/lib.bash"
ECOLI_GENOMES=$B/data/ecoli/genomes
PNEUMO_GENOMES=$B/data/pneumo/genomes
ECOLI_FULL=${ECOLI_FULL:-$B/indexes/ecoli/scale_3682.thm2}
PNEUMO_FULL=$B/indexes/pneumo/s_pneumoniae.thm2
read -r -a DISJOINT <<< "${DISJOINT:-1 36 368 1841}"
read -r -a SHARED <<< "${SHARED:-10 50 100}"
read -r -a PARTS <<< "${PARTS:-2 4 8 16}"
PNEUMO_SPLIT=(121033 60517)
W_ROOT=${W_ROOT:-$B/work}
IDX_ROOT=${IDX_ROOT:-$B/indexes}
FOF=$B/work/fof
SETS=${SETS:-$P/sets/ecoli}
mkdir -p "$FOF" "$IDX_ROOT/ecoli/merge" "$IDX_ROOT/pneumo/merge"
[[ -f $ECOLI_FULL ]] || { log "no $ECOLI_FULL"; exit 1; }

# direct <prefix> <fof> <output> <round> <cap>: the direct build, as in queue_subset.bash
direct() {
    local pre=$1 f=$2 out=$3 r=$4 g=$W/$1.genomes.txt
    local CAP_GB=$5
    done_run "${pre}_themisto_r$r" && { log "skip $pre round $r (in the ledger)"; return 0; }
    [[ -s $f ]] || { log "FAILED direct build $pre: empty genome list $f"; return 1; }
    cut -f1 "$f" > "$g"
    INPUTS="@$f" run "${pre}_ggcat_r$r" "$W/$pre.unitigs.fa" 0 \
        "$BIN/ggcat" build -l "$g" -s 1 -k $K -m 8 -p -j "$THREADS" -t "$W/tmp" -o "$W/$pre.unitigs.fa" &&
    INPUTS="$W/$pre.unitigs.fa" run "${pre}_sbwt_r$r" "$W/$pre.sbwt" 0 \
        "$BIN/sbwt" build -i "$W/$pre.unitigs.fa" -o "$W/$pre" -k $K -t "$THREADS" -m 8 -r -v --temp-dir "$W/tmp" &&
    INPUTS="@$f $W/$pre.sbwt" run "${pre}_themisto_r$r" "$out" 0 \
        "$BIN/themisto2" build --file-colors "$f" -s "$W/$pre.sbwt" -o "$out" -k $K -d 30 -t "$THREADS" --temp-dir "$W/tmp" --no-lcs ||
        { log "FAILED direct build $pre round $r"; return 1; }
    rm -f -- "$W/$pre.unitigs.fa" "$W/$pre.sbwt" "$g"
}
# compare <description> <index1> <index2>: check, once (skipped if in checks.tsv), when both exist
compare() {
    [[ -f $CHECKS ]] && awk -F'\t' -v d="$1" '$1 == d { f = 1 } END { exit !f }' "$CHECKS" && return 0
    [[ -f $2 && -f $3 ]] || { log "no check of '$1': $2 or $3 missing"; return 1; }
    check "$CHECKS" "$1" "$2" "$3"
}
# compare_on_subset <description> <index> <full index> <colors file> <full index cut to the colors>:
# check of <index> against <full index> on the given colors only (both cut by subset, capped, not timed;
# the cut of the full index is made once and kept). Skipped if in checks.tsv.
compare_on_subset() {
    local d=$1 m=$2 full=$3 colors=$4 ref=$5 cut=$W/check_cut.thm2
    local capped=(systemd-run --user --scope --quiet -p MemoryMax="${SUBSET_CHECK_CAP:-160}G" -p MemorySwapMax=0 nice -n 19)
    [[ -f $CHECKS ]] && awk -F'\t' -v d="$d" '$1 == d { f = 1 } END { exit !f }' "$CHECKS" && return 0
    [[ -f $m ]] || { log "no check of '$d': $m missing"; return 1; }
    if [[ ! -f $ref ]]; then
        "${capped[@]}" "$BIN/themisto2" subset --index "$full" --colors-file "$colors" -o "${ref%.thm2}.part.thm2" -t "$THREADS" --no-lcs 2> "$W/check_cut_full.log" &&
            mv "${ref%.thm2}.part.thm2" "$ref" || { log "cut of the full index failed: $d"; rm -f -- "${ref%.thm2}.part.thm2"; return 1; }
    fi
    "${capped[@]}" "$BIN/themisto2" subset --index "$m" --colors-file "$colors" -o "$cut" -t "$THREADS" --no-lcs 2> "$W/check_cut.log" ||
        { log "cut of the result failed: $d"; rm -f -- "$cut"; return 1; }
    CAP_GB=16 check "$CHECKS" "$d" "$cut" "$ref"
    rm -f -- "$cut"
}
# input <name> <full index> <colors file> <output> <cap>: an input index by subset (made once, kept)
input() {
    [[ -f $4 ]] && return 0
    CAP_GB=$5 INPUTS="$2" run "$1" "$4" 0 "$BIN/themisto2" subset --index "$2" --colors-file "$3" -o "$4" -t "$THREADS" --no-lcs ||
        { log "FAILED input $1"; return 1; }
}
# merge_run <name> <round> <cap> <output> <flags> <index>...: one timed merge (its output is deleted by
# the caller: the one-by-one merges need it for their next step)
merge_run() {
    local name=$1 r=$2 cap=$3 out=$4 flags=$5; shift 5
    local list=$W/$name.list
    printf '%s\n' "$@" > "$list"
    CAP_GB=$cap INPUTS="$*" run "${name}_r$r" "$out" 0 \
        "$BIN/themisto2" merge --index-file-list "$list" --temp-dir "$W/tmp" -o "$out" -t "$THREADS" --no-lcs $flags ||
        { log "FAILED ${name}_r$r"; return 1; }
}
# merged <dataset> <name> <round> <cap> <check cap> <reference> <flags> <index>...: merge_run, then in
# round 1 the check against the reference, then the result is deleted
merged() {
    local ds=$1 name=$2 r=$3 cap=$4 ccap=$5 ref=$6 flags=$7; shift 7
    local out=$W/$name.thm2
    merge_run "${ds}_merge_$name" "$r" "$cap" "$out" "$flags" "$@"
    if (( r == 1 )); then
        if [[ $ds == pneumo ]]; then
            compare_on_subset "$ds merge $name = full index, on keep_1000" "$out" "$ref" "$P/sets/pneumo/keep_1000.txt" "$PI/full_keep_1000.thm2"
        else
            CAP_GB=$ccap compare "$ds merge $name = full index" "$out" "$ref"
        fi
    fi
    rm -f -- "$out"
}
# plain <dataset> <set> <fof> <round>: plain themisto2 build from the genomes (small sets of queue_subset.bash)
plain() {
    local out=$W/plain_$2.thm2
    CAP_GB=16 INPUTS="@$3" run "${1}_plainbuild_$2_r$4" "$out" 1 \
        "$BIN/themisto2" build --file-colors "$3" -o "$out" -k $K -d 30 -t "$THREADS" --temp-dir "$W/tmp" --no-lcs || log "FAILED ${1}_plainbuild_$2_r$4"
}

# ---- E. coli
ecoli_fof() { awk -v g="$ECOLI_GENOMES" '{ print g "/" $1 ".fna\t" $1 }' "$2" > "$FOF/ecoli_$1.tsv"; echo "$FOF/ecoli_$1.tsv"; }
EI=$IDX_ROOT/ecoli/merge
ecoli_inputs() {
    local m p k i
    W=$W_ROOT/ecoli_merge; mkdir -p "$W/tmp"
    for m in "${DISJOINT[@]}"; do for p in A B; do
        input "ecoli_input_disjoint_${m}_$p" "$ECOLI_FULL" "$SETS/merge_disjoint_${m}_$p.txt" "$EI/disjoint_${m}_$p.thm2" 16 || return 1
    done; done
    for m in "${SHARED[@]}"; do for p in A B; do
        input "ecoli_input_shared_${m}_$p" "$ECOLI_FULL" "$SETS/merge_shared_${m}_$p.txt" "$EI/shared_${m}_$p.thm2" 16 || return 1
    done; done
    for k in "${PARTS[@]}"; do for (( i = 0; i < k; i++ )); do
        input "ecoli_input_parts${k}_$i" "$ECOLI_FULL" "$SETS/parts${k}_$i.txt" "$EI/parts${k}_$i.thm2" 16 || return 1
    done; done
}
# one_by_one <K> <round>: ((p0 + p1) + p2) + ...; the intermediates are kept until the last step is done
one_by_one() {
    local k=$1 r=$2 i prev=$EI/parts$1_0.thm2 out name=ecoli_merge_onebyone$1
    if ! done_run "${name}_s$(( k - 1 ))_r$r"; then
        for (( i = 1; i < k; i++ )); do
            out=$W/onebyone$k.$i.thm2
            merge_run "${name}_s$i" "$r" 16 "$out" "" "$prev" "$EI/parts${k}_$i.thm2" || return 1
            prev=$out
        done
    fi
    (( r == 1 )) && CAP_GB=16 compare "ecoli merge onebyone$k = full index" "$W/onebyone$k.$(( k - 1 )).thm2" "$ECOLI_FULL"
    rm -f -- "$W"/onebyone$k.*.thm2
}
ecoli_round() {
    local r=$1 m k i
    W=$W_ROOT/ecoli_merge; mkdir -p "$W/tmp"
    for m in "${DISJOINT[@]}"; do
        merged ecoli "disjoint_$m" "$r" 16 16 "$ECOLI_FULL" "" "$EI/disjoint_${m}_A.thm2" "$EI/disjoint_${m}_B.thm2"
    done
    for m in "${SHARED[@]}"; do
        merged ecoli "shared_$m" "$r" 16 16 "$ECOLI_FULL" "--merge-shared-colors" "$EI/shared_${m}_A.thm2" "$EI/shared_${m}_B.thm2"
    done
    for k in "${PARTS[@]}"; do
        local parts=(); for (( i = 0; i < k; i++ )); do parts+=("$EI/parts${k}_$i.thm2"); done
        merged ecoli "multi$k" "$r" 16 16 "$ECOLI_FULL" "" "${parts[@]}"
        (( k > 2 )) && one_by_one "$k" "$r"
    done
    m=${DISJOINT[-1]}
    merged ecoli "disjoint_${m}_keepdummies" "$r" 16 16 "$ECOLI_FULL" "--keep-redundant-dummies" "$EI/disjoint_${m}_A.thm2" "$EI/disjoint_${m}_B.thm2"
    merged ecoli "disjoint_${m}_lowram" "$r" 16 16 "$ECOLI_FULL" "--low-ram-mode" "$EI/disjoint_${m}_A.thm2" "$EI/disjoint_${m}_B.thm2"
    # The direct build of the union (round 1: the build of the input index, by queue_subset.bash)
    (( r > 1 )) && { direct ecoli_full "$(ecoli_fof all "$SETS/order.txt")" "$W/full_rebuild.thm2" "$r" 16; rm -f -- "$W/full_rebuild.thm2"; }
    for k in 1 10 100; do plain ecoli "keep_$k" "$(ecoli_fof "keep_$k" "$SETS/keep_$k.txt")" "$r"; done
}

# ---- S. pneumoniae
PI=$IDX_ROOT/pneumo/merge
NAMES=$P/sets/pneumo/color_names.tsv
pneumo_inputs() {
    local n
    W=$W_ROOT/pneumo_merge; mkdir -p "$W/tmp"
    for n in "${PNEUMO_SPLIT[@]}"; do
        cut -f2 "$NAMES" | head -n "$n" > "$FOF/pneumo_first_$n.txt"
        cut -f2 "$NAMES" | tail -n +$(( n + 1 )) > "$FOF/pneumo_rest_$n.txt"
        input "pneumo_input_first_$n" "$PNEUMO_FULL" "$FOF/pneumo_first_$n.txt" "$PI/first_$n.thm2" 160 || return 1
        input "pneumo_input_rest_$n" "$PNEUMO_FULL" "$FOF/pneumo_rest_$n.txt" "$PI/rest_$n.thm2" 160 || return 1
    done
}
pneumo_fof() {
    awk -F'\t' -v g="$PNEUMO_GENOMES" 'NR == FNR { s[$1] = 1; next } ($2 in s) { a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' \
        "$P/sets/pneumo/$1.txt" "$NAMES" > "$FOF/pneumo_$1.tsv"
    echo "$FOF/pneumo_$1.tsv"
}
pneumo_round() {
    local r=$1 n k
    W=$W_ROOT/pneumo_merge; mkdir -p "$W/tmp"
    for n in "${PNEUMO_SPLIT[@]}"; do
        merged pneumo "first_$n" "$r" "$MERGE_CAP_PNEUMO" 160 "$PNEUMO_FULL" "" "$PI/first_$n.thm2" "$PI/rest_$n.thm2"
    done
    for k in 1 10 100; do plain pneumo "keep_$k" "$(pneumo_fof "keep_$k")" "$r"; done
}

ecoli_inputs || exit 1
(( PNEUMO )) && { pneumo_inputs || exit 1; }
for (( r = 1; r <= REPS_ECOLI || (PNEUMO && r <= REPS_PNEUMO); r++ )); do
    log "round $r"
    (( r <= REPS_ECOLI )) && ecoli_round "$r"
    (( PNEUMO && r <= REPS_PNEUMO )) && pneumo_round "$r"
done
log "queue_merge done"
