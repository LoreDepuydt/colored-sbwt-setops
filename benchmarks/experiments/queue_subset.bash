#!/bin/bash
# The subset benchmark: subset against building the same index directly, on E. coli and
# S. pneumoniae. Every run goes through run.bash (lock, quiet server, cap, cold start, samples,
# ledger) and is skipped if it is in the ledger already, so the queue can be stopped and restarted.
#
# Repetitions round by round: round r runs every E. coli run once, then (r <= REPS_PNEUMO) every
# pneumo run once, so that a first complete pass is done early. Run names end in _r<round>.
#
# E. coli (input: the full index of the 3,682 genomes, built once by the direct pipeline below):
#   subset keep_<n> and its complement (--complement), n = 1, 10, 100, 368 (10 %), 1841 (50 %), 3314 (90 %)
#   the direct build of the same genomes, for each of these 12 sets
#   keep_1841 (SCALING_N) also with 1, 2, 4, 8, 16 threads (thread scaling) and with --low-ram-mode
# S. pneumoniae (input: s_pneumoniae.thm2, 121,034 genomes):
#   subset keep_<n> and its complement, n = 1, 10, 100, 1000, 12103 (10 %), 60517 (50 %), 108931 (90 %)
#   the direct build of keep_<n> for n <= 12103 (the ATB genomes downloaded by download_pneumo.bash)
# Caps: E. coli 16 GiB; pneumo subsets 160 GiB (at 100 GiB the complements, 80 GiB of anonymous memory,
# hit the cap with the file cache), pneumo builds 64 GiB.
# Direct build: the fastest pipeline that we found: ggcat build
# (no colors, -s 1) -> sbwt build -r -> themisto2 build --sbwt --file-colors; three runs
# <prefix>_ggcat/_sbwt/_themisto. The genomes are given in the order of the colors in the input
# index, so that the subset and the build number their colors the same way.
# Correctness, in round 1: each subset against its direct build (check in lib.bash: the k-mers, the
# color set of every k-mer, and the color names), one line in checks.tsv.
#
# Usage: bash queue_subset.bash   (settings: THREADS (32), REPS_ECOLI (5), REPS_PNEUMO (3); for a small
#   test: ECOLI_FULL_SET (a set file, default sets/ecoli/order.txt), ECOLI_FULL, ECOLI_N, SCALING_N, W_ROOT)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS_ECOLI=${REPS_ECOLI:-5}
REPS_PNEUMO=${REPS_PNEUMO:-3}
CHECKS=${CHECKS:-$P/checks.tsv}
source "$P/lib.bash"
ECOLI_GENOMES=$B/data/ecoli/genomes
PNEUMO_GENOMES=$B/data/pneumo/genomes
ECOLI_FULL_SET=${ECOLI_FULL_SET:-$P/sets/ecoli/order.txt}
ECOLI_FULL=${ECOLI_FULL:-$B/indexes/ecoli/scale_3682.thm2}
PNEUMO_FULL=$B/indexes/pneumo/s_pneumoniae.thm2
read -r -a ECOLI_N <<< "${ECOLI_N:-1 10 100 368 1841 3314}"
SCALING_N=${SCALING_N:-1841}
W_ROOT=${W_ROOT:-$B/work}
PNEUMO_N=(1 10 100 1000 12103 60517 108931)
PNEUMO_BUILD_MAX=12103
FOF=$B/work/fof
[[ -s $ECOLI_FULL_SET ]] || { log "no set $ECOLI_FULL_SET"; exit 1; }
mkdir -p "$(dirname "$ECOLI_FULL")" "$FOF"

# direct <prefix> <fof: genome TAB color name> <output index> <round> <cap GiB>: the direct build;
# its intermediates are deleted only once the build is done, so that a restart can go on from them
direct() {
    local pre=$1 f=$2 out=$3 r=$4 g=$W/$1.genomes.txt
    local CAP_GB=$5
    # Done already: return before the cleanup below, which would delete the intermediates of the
    # same set in another round (they have the same names)
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
# one_set <dataset> <set name> <subset flags> <full index> <fof of the direct build, or ""> <round> <subset cap> <build cap>
one_set() {
    local ds=$1 s=$2 flags=$3 full=$4 f=$5 r=$6 sub=$W/$2.subset.thm2 dir=$W/$2.direct.thm2
    # In round 1 the results are kept until the check; later rounds delete them after each run
    local del=$(( r > 1 )) CAP_GB=$7
    INPUTS="$full" run "${ds}_subset_${s}_r$r" "$sub" "$del" \
        "$BIN/themisto2" subset --index "$full" $flags -o "$sub" -t "$THREADS" --no-lcs || log "FAILED ${ds}_subset_${s}_r$r"
    if [[ -n $f ]]; then
        direct "${ds}_direct_$s" "$f" "$dir" "$r" "$8"
        (( r == 1 )) && CAP_GB=$8 compare "$ds $s: subset = direct build" "$sub" "$dir"
        rm -f -- "$dir"
    fi
    rm -f -- "$sub"
}

# ---- E. coli
ecoli_fof() { awk -v g="$ECOLI_GENOMES" '{ print g "/" $1 ".fna\t" $1 }' "$2" > "$FOF/ecoli_$1.tsv"; echo "$FOF/ecoli_$1.tsv"; }
ecoli_round() {
    local r=$1 n
    W=$W_ROOT/ecoli; mkdir -p "$W/tmp"
    for n in "${ECOLI_N[@]}"; do
        # keep_<n>: the first n of the random order; compl_<n>: the others (the full index has the order of order.txt)
        tail -n +$(( n + 1 )) "$ECOLI_FULL_SET" > "$FOF/ecoli_compl_$n.txt"
        one_set ecoli "keep_$n" "--colors-file $P/sets/ecoli/keep_$n.txt" "$ECOLI_FULL" "$(ecoli_fof "keep_$n" "$P/sets/ecoli/keep_$n.txt")" "$r" 16 16
        one_set ecoli "compl_$n" "--complement --colors-file $P/sets/ecoli/keep_$n.txt" "$ECOLI_FULL" "$(ecoli_fof "compl_$n" "$FOF/ecoli_compl_$n.txt")" "$r" 16 16
    done
    local t s=$P/sets/ecoli/keep_$SCALING_N.txt o=$W/scaling.thm2
    for t in 1 2 4 8 16; do
        CAP_GB=16 THREADS=$t INPUTS="$ECOLI_FULL" run "ecoli_subset_keep_${SCALING_N}_t${t}_r$r" "$o" 1 \
            "$BIN/themisto2" subset --index "$ECOLI_FULL" --colors-file "$s" -o "$o" -t "$t" --no-lcs || log "FAILED t$t round $r"
    done
    CAP_GB=16 INPUTS="$ECOLI_FULL" run "ecoli_subset_keep_${SCALING_N}_lowram_r$r" "$o" 1 \
        "$BIN/themisto2" subset --index "$ECOLI_FULL" --colors-file "$s" -o "$o" -t "$THREADS" --no-lcs --low-ram-mode || log "FAILED lowram round $r"
}

# ---- S. pneumoniae
# The genomes of a set with the color names of the index, in the order of the index
pneumo_fof() {
    awk -F'\t' -v g="$PNEUMO_GENOMES" 'NR == FNR { s[$1] = 1; next } ($2 in s) { a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' \
        "$P/sets/pneumo/$1.txt" "$P/sets/pneumo/color_names.tsv" > "$FOF/pneumo_$1.tsv"
    echo "$FOF/pneumo_$1.tsv"
}
pneumo_round() {
    local r=$1 n f
    W=$W_ROOT/pneumo; mkdir -p "$W/tmp"
    for n in "${PNEUMO_N[@]}"; do
        f=""; (( n <= PNEUMO_BUILD_MAX )) && f=$(pneumo_fof "keep_$n")
        one_set pneumo "keep_$n" "--colors-file $P/sets/pneumo/keep_$n.txt" "$PNEUMO_FULL" "$f" "$r" 160 64
        one_set pneumo "compl_$n" "--complement --colors-file $P/sets/pneumo/keep_$n.txt" "$PNEUMO_FULL" "" "$r" 160 64
    done
}

# ---- The input: the full E. coli index (one direct build, kept)
W=$W_ROOT/ecoli; mkdir -p "$W/tmp"
if [[ ! -f $ECOLI_FULL ]]; then
    direct ecoli_full "$(ecoli_fof all "$ECOLI_FULL_SET")" "$ECOLI_FULL" 1 16 || exit 1
fi
[[ -f $PNEUMO_FULL ]] || { log "no $PNEUMO_FULL"; exit 1; }
for f in $(pneumo_fof "keep_$PNEUMO_BUILD_MAX" | xargs cut -f1); do [[ -f $f ]] || { log "missing $f: bash download_pneumo.bash keep_$PNEUMO_BUILD_MAX"; exit 1; }; done

for (( r = 1; r <= REPS_ECOLI || r <= REPS_PNEUMO; r++ )); do
    log "round $r"
    (( r <= REPS_ECOLI )) && ecoli_round "$r"
    (( r <= REPS_PNEUMO )) && pneumo_round "$r"
done
log "queue_subset done"
