#!/bin/bash
# The build time of the parts that are actually merged (the merge rows of the comparison with a
# rebuild), instead of random sets of the same size. Every run through run.bash,
# skipped if it is in the ledger, round by round:
#   E. coli: sets/ecoli/merge_disjoint_1_B.txt (1 genome) and merge_disjoint_368_B.txt (368), REPS_ECOLI (5)
#   S. pneumoniae: the 1 genome that "merge: add 1 genome" adds (rest_121033: the last color of
#     s_pneumoniae.thm2, SAMN42364221), REPS_PNEUMO (3)
# Each part: the rebuild pipeline of the other rows (direct in queue_merge.bash: ggcat -> sbwt -> themisto2
# build --sbwt, runs <ds>_partdirect_<part>_{ggcat,sbwt,themisto}_r<i>), and for at most 100 genomes also
# plain themisto2 build (runs <ds>_partplain_<part>_r<i>). Caps 16 GiB, as the builds of queue_merge.bash.
# The outputs are deleted after each run. Summary: python3 parts_summary.py.
# Usage: bash queue_parts.bash   (settings: THREADS (32), REPS_ECOLI, REPS_PNEUMO; for a test: LEDGER,
#   RESULTS, PARTS ("ecoli:merge_disjoint_1_B pneumo:rest_121033"))
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
REPS_ECOLI=${REPS_ECOLI:-5}
REPS_PNEUMO=${REPS_PNEUMO:-3}
source "$P/lib.bash"
PARTS=${PARTS:-"ecoli:merge_disjoint_1_B ecoli:merge_disjoint_368_B pneumo:rest_121033"}
W=$B/work/parts; mkdir -p "$W/tmp"
CAP_GB=16
NAMES=$P/sets/pneumo/color_names.tsv

# fof <dataset> <part>: genome TAB color name, as queue_merge.bash
fof() {
    local f=$W/$1_$2.fof.tsv
    if [[ $1 == ecoli ]]; then
        awk -v g="$B/data/ecoli/genomes" '{ print g "/" $1 ".fna\t" $1 }' "$P/sets/ecoli/$2.txt" > "$f"
    else  # rest_<n>: the colors of s_pneumoniae.thm2 after the first n
        local n=${2#rest_}
        tail -n +$(( n + 1 )) "$NAMES" | awk -F'\t' -v g="$B/data/pneumo/genomes" '{ a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' > "$f"
    fi
    echo "$f"
}
direct() {
    local pre=$1 f=$2 r=$3 g=$W/$1.genomes.txt out=$W/$1.thm2
    done_run "${pre}_themisto_r$r" && { log "skip $pre round $r (in the ledger)"; return 0; }
    cut -f1 "$f" > "$g"
    INPUTS="@$f" run "${pre}_ggcat_r$r" "$W/$pre.unitigs.fa" 0 \
        "$BIN/ggcat" build -l "$g" -s 1 -k $K -m 8 -p -j "$THREADS" -t "$W/tmp" -o "$W/$pre.unitigs.fa" &&
    INPUTS="$W/$pre.unitigs.fa" run "${pre}_sbwt_r$r" "$W/$pre.sbwt" 0 \
        "$BIN/sbwt" build -i "$W/$pre.unitigs.fa" -o "$W/$pre" -k $K -t "$THREADS" -m 8 -r -v --temp-dir "$W/tmp" &&
    INPUTS="@$f $W/$pre.sbwt" run "${pre}_themisto_r$r" "$out" 1 \
        "$BIN/themisto2" build --file-colors "$f" -s "$W/$pre.sbwt" -o "$out" -k $K -d 30 -t "$THREADS" --temp-dir "$W/tmp" --no-lcs ||
        log "FAILED direct build $pre round $r"
    rm -f -- "$W/$pre.unitigs.fa" "$W/$pre.sbwt" "$g" "$out"
}
plain() {
    local out=$W/$1.plain.thm2
    INPUTS="@$2" run "$1_r$3" "$out" 1 \
        "$BIN/themisto2" build --file-colors "$2" -o "$out" -k $K -d 30 -t "$THREADS" --temp-dir "$W/tmp" --no-lcs || log "FAILED $1_r$3"
    rm -f -- "$out"
}

for p in $PARTS; do
    ds=${p%%:*} part=${p#*:}; f=$(fof "$ds" "$part")
    [[ -s $f ]] || { log "ABORT: empty genome list for $p"; exit 1; }
    while IFS=$'\t' read -r g _; do [[ -s $g ]] || { log "ABORT: missing genome $g"; exit 1; }; done < "$f"
done
max=$(( REPS_ECOLI > REPS_PNEUMO ? REPS_ECOLI : REPS_PNEUMO ))
for (( r = 1; r <= max; r++ )); do
    log "round $r"
    for p in $PARTS; do
        ds=${p%%:*} part=${p#*:}
        [[ $ds == ecoli ]] && reps=$REPS_ECOLI || reps=$REPS_PNEUMO
        (( r <= reps )) || continue
        f=$W/${ds}_$part.fof.tsv
        direct "${ds}_partdirect_$part" "$f" "$r"
        (( $(wc -l < "$f") <= 100 )) && plain "${ds}_partplain_$part" "$f" "$r"
    done
done
log "queue_parts done"
