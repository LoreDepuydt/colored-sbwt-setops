#!/bin/bash
# The full S. pneumoniae rebuilds, the only rebuilds the subset and merge of
# one genome had not: the direct build (the pipeline of queue_subset.bash: ggcat -> sbwt -> themisto2
# build --sbwt, k = 31, --no-lcs) of
#   compl_1  all genomes but the first of sets/pneumo/order.txt (121,033): rebuild for "subset: remove 1"
#            runs pneumo_direct_compl_1_{ggcat,sbwt,themisto}_r<i>
#   all      all 121,034 genomes: rebuild for "merge: add 1"; runs pneumo_full_{ggcat,sbwt,themisto}_r<i>
# Every run through run.bash, skipped if in the ledger. REPS repetitions (default 1). Cap CAP GiB (200)
# for every step. The genomes (download_pneumo.bash all) are given in the order of the colors of
# s_pneumoniae.thm2, so that the rebuild numbers its colors the same way.
# Checks (round 1, into checks.tsv; under the lock of run.bash), without a full export (829 GB for this index):
#   - compl_1: the rebuild has as many k-mers as pneumo_subset_compl_1 kept (414,072,812)
#   - all: themisto2 stats of the rebuild = of s_pneumoniae.thm2 (the fields of STATS_FIELDS; both
#     stats are kept in the run folder of the themisto step), and the rebuild cut to keep_1000 = the
#     index cut to keep_1000 (indexes/pneumo/merge/full_keep_1000.thm2), by check of lib.bash
# The rebuild of all is kept ($W/all.direct.thm2) for a full comparison later; compl_1 is deleted after
# its check. A later round overwrites both.
# Usage: bash queue_rebuild.bash   (settings: THREADS (32), REPS (1), CAP (200); for a small test: BUILDS,
#   COMPL_PRE, COMPL_KMERS, FULL_PRE, REF, REF_CUT ("" = cut REF here), CHECK_SET, W_ROOT)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS=${REPS:-1}
CAP=${CAP:-200}
CHECKS=${CHECKS:-$P/checks.tsv}
source "$P/lib.bash"
PNEUMO_GENOMES=$B/data/pneumo/genomes
W_ROOT=${W_ROOT:-$B/work}
W=$W_ROOT/pneumo_rebuild; mkdir -p "$W/tmp"
# <set>:<run prefix>
BUILDS=${BUILDS:-"compl_1:pneumo_direct_compl_1 all:pneumo_full"}
COMPL_PRE=${COMPL_PRE:-pneumo_direct_compl_1}
COMPL_KMERS=${COMPL_KMERS:-414072812}
FULL_PRE=${FULL_PRE:-pneumo_full}
REF=${REF:-$B/indexes/pneumo/s_pneumoniae.thm2}
CHECK_SET=${CHECK_SET:-keep_1000}
REF_CUT=${REF_CUT-$B/indexes/pneumo/merge/full_keep_1000.thm2}
STATS_FIELDS="n_kmers n_colors sbwt_length n_sampled_positions n_distinct_color_sets n_sparse_color_sets n_dense_color_sets color_set_storage_size_bytes n_unitigs"

# direct <prefix> <fof: genome TAB color name> <output index> <round> <cap GiB>: as in queue_subset.bash
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
# The genomes of a set with the color names of the index, in the order of the index (as queue_subset.bash)
fof() {
    awk -F'\t' -v g="$PNEUMO_GENOMES" 'NR == FNR { s[$1] = 1; next } ($2 in s) { a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' \
        "$P/sets/pneumo/$1.txt" "$P/sets/pneumo/color_names.tsv" > "$W/$1.fof.tsv"
    echo "$W/$1.fof.tsv"
}
capped() { systemd-run --user --scope --quiet -p MemoryMax="${CAP}G" -p MemorySwapMax=0 nice -n 19 "$@"; }
done_check() { [[ -f $CHECKS ]] && awk -F'\t' -v d="$1" '$1 == d { f = 1 } END { exit !f }' "$CHECKS"; }
# record <description> <result>: one line in checks.tsv (no color names column)
record() {
    [[ -s $CHECKS ]] || printf 'check\tkmers_and_color_sets\tcolor_names\tdate\n' > "$CHECKS"
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "n/a" "$(date -Is)" | tee -a "$CHECKS"
}
field() { awk -F'\t' -v k="$1" '$1 == k { print $2 }' "$2"; }

check_compl() {
    local out=$1 run=$2 d="pneumo compl_1: the rebuild has the k-mers kept by the subset (n_kmers = $COMPL_KMERS)" n
    done_check "$d" && return 0
    [[ -f $out ]] || { log "no check of $out: missing"; return 1; }
    capped "$BIN/themisto2" stats -i "$out" -t "$THREADS" 2> /dev/null > "$run/stats.tsv" || { log "no check of $out: stats failed"; return 1; }
    n=$(field n_kmers "$run/stats.tsv")
    record "$d" "$( [[ $n == "$COMPL_KMERS" ]] && echo equal || echo "DIFFERENT ($n)")"
}
check_full() {
    local out=$1 run=$2 d1="pneumo all: themisto2 stats of the rebuild = of $(basename "$REF")" diff="" k a b
    local d2="pneumo all: rebuild = $(basename "$REF"), both cut to $CHECK_SET" cut=$W/rebuild_$CHECK_SET.thm2 ref_cut=$REF_CUT
    done_check "$d1" && done_check "$d2" && return 0
    [[ -f $out ]] || { log "no check of $out: missing"; return 1; }
    if ! done_check "$d1"; then
        capped "$BIN/themisto2" stats -i "$out" -t "$THREADS" 2> /dev/null > "$run/stats.tsv" &&
            capped "$BIN/themisto2" stats -i "$REF" -t "$THREADS" 2> /dev/null > "$run/stats_ref.tsv" ||
            { log "no check of $out: stats failed"; return 1; }
        for k in $STATS_FIELDS; do
            a=$(field "$k" "$run/stats.tsv"); b=$(field "$k" "$run/stats_ref.tsv")
            [[ $a == "$b" ]] || diff+="${diff:+; }$k $a vs $b"
        done
        record "$d1" "$( [[ -z $diff ]] && echo equal || echo "DIFFERENT ($diff)")"
    fi
    done_check "$d2" && return 0
    capped "$BIN/themisto2" subset --index "$out" --colors-file "$P/sets/pneumo/$CHECK_SET.txt" -o "$cut" -t "$THREADS" --no-lcs 2> "$W/cut.log" ||
        { log "no check of $out: cutting it to $CHECK_SET failed (see $W/cut.log)"; return 1; }
    if [[ -z $ref_cut ]]; then
        ref_cut=$W/ref_$CHECK_SET.thm2
        capped "$BIN/themisto2" subset --index "$REF" --colors-file "$P/sets/pneumo/$CHECK_SET.txt" -o "$ref_cut" -t "$THREADS" --no-lcs 2> "$W/cut_ref.log" ||
            { log "no check of $out: cutting $REF to $CHECK_SET failed"; return 1; }
    fi
    CAP_GB=$CAP check "$CHECKS" "$d2" "$cut" "$ref_cut"
    rm -f -- "$cut" "$W/ref_$CHECK_SET.thm2"
}

for b in $BUILDS; do
    s=${b%%:*}
    for f in $(fof "$s" | xargs cut -f1); do [[ -f $f ]] || { log "missing $f: bash download_pneumo.bash $s"; exit 1; }; done
done
for (( r = 1; r <= REPS; r++ )); do
    log "round $r"
    for b in $BUILDS; do
        s=${b%%:*} pre=${b#*:} out=$W/$s.direct.thm2
        done_run "${pre}_themisto_r$r" || rm -f -- "$out"
        direct "$pre" "$(fof "$s")" "$out" "$r" "$CAP" || continue
        (( r == 1 )) || { [[ $pre == "$COMPL_PRE" ]] && rm -f -- "$out"; continue; }
        run_dir=${RESULTS:-$P/runs}/${pre}_themisto_r1_$(hostname -s)
        # The checks load whole indexes (stats, a cut of the full index): under the lock of run.bash too
        if [[ $pre == "$COMPL_PRE" ]]; then
            ( flock 9; check_compl "$out" "$run_dir" ) 9> "$P/.lock" && rm -f -- "$out"
        elif [[ $pre == "$FULL_PRE" ]]; then
            ( flock 9; check_full "$out" "$run_dir" ) 9> "$P/.lock"
        fi
    done
done
log "queue_rebuild done"
