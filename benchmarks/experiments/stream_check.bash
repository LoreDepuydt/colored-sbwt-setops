#!/bin/bash
# Full checks of S. pneumoniae results (whole index, not only keep_1000), with the comparison of check in
# lib.bash: compare_color_dumps (the k-mers and the color set of every k-mer, by color id) and the color
# names. A whole S. pneumoniae export is over 800 GB, so the color sets are never written uncompressed:
# themisto2 export --colors-to-stdout | zstd, and compare_color_dumps reads them back through a FIFO
# (it reads each color_sets.txt once, front to back, after both metadata and unitig files). The
# export of the rebuild is made once and used for both merges.
#
# Here:
#   - the merges "add 1 genome" (first_121033) and "two halves" (first_60517) are run again (timed,
#     through run.bash, as round 4, inputs indexes/pneumo/merge/{first,rest}_<n>.thm2), their results
#     are kept in work/pneumo_fullcheck/, and each is compared with the full rebuild
#     work/pneumo_rebuild/all.direct.thm2
#   - the subset compl_1 against its rebuild pneumo_direct_compl_1: that rebuild was deleted after its
#     n_kmers check, so it is not done here
# Everything that loads an index (the merges through run.bash, exports, compares) is under the shared lock.
# Disk: stops if the disk of $B falls under CHECK_MIN_FREE_GB (300) GiB free. Each compressed export is
# about 1/8 of the text (keep_1000: 1.24 GB -> 157 MB), so about 150 GB for a whole index. The export of a
# merge result is deleted (by name) after its compare, the export of the rebuild after both; the merge
# results too, unless KEEP_MERGES=1.
# Usage: bash stream_check.bash   (settings: THREADS (32), CAP (240 for the merges), CHECK_CAP (100);
#   for a small test: PAIRS ("<description>|<index A>|<index B>" ...), MERGES=0, W_ROOT, CHECKS)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
CAP=${CAP:-240}
CHECK_CAP=${CHECK_CAP:-100}
MERGES=${MERGES:-1}
KEEP_MERGES=${KEEP_MERGES:-0}
ROUND=${ROUND:-4}
CHECKS=${CHECKS:-$P/checks.tsv}
source "$P/lib.bash"
W_ROOT=${W_ROOT:-$B/work}
W=$W_ROOT/pneumo_fullcheck; X=$W/export; mkdir -p "$W/tmp" "$X"
PI=$B/indexes/pneumo/merge
REBUILD=$B/work/pneumo_rebuild/all.direct.thm2
LOCK=$P/.lock
MIN_FREE=${CHECK_MIN_FREE_GB:-300}
capped() { systemd-run --user --scope --quiet -p MemoryMax="${CHECK_CAP}G" -p MemorySwapMax=0 nice -n 19 "$@"; }
free_gb() { df --output=avail -B1G "$B" | tail -1; }
done_check() { [[ -f $CHECKS ]] && awk -F'\t' -v d="$1" '$1 == d { f = 1 } END { exit !f }' "$CHECKS"; }

# export_z <index> <prefix>: <prefix>.metadata.txt, <prefix>.unitigs.fa and <prefix>.color_sets.txt.zst
# (made once: skipped if <prefix>.done exists). Stopped (whole process group) if the disk runs low.
export_z() {
    local i=$1 x=$2 pid st
    [[ -f $x.done ]] && { log "export of $i: kept from before ($x)"; return 0; }
    (( $(free_gb) >= MIN_FREE )) || { log "no export of $i: under $MIN_FREE GiB free"; return 1; }
    log "export $i -> $x"
    setsid bash -c 'set -o pipefail; systemd-run --user --scope --quiet -p MemoryMax="$1G" -p MemorySwapMax=0 nice -n 19 \
        "$2" export -i "$3" -o "$4" -t "$5" --colors-to-stdout 2> "$4.log" | nice -n 19 zstd -q -T8 -3 -o "$4.color_sets.txt.zst.part" -f' \
        _ "$CHECK_CAP" "$BIN/themisto2" "$i" "$x" "$THREADS" & pid=$!
    while kill -0 "$pid" 2> /dev/null; do
        if (( $(free_gb) < MIN_FREE )); then
            kill -TERM -- "-$pid"; wait "$pid"; rm -f -- "$x.color_sets.txt.zst.part"; log "export of $i stopped: under $MIN_FREE GiB free"; return 1
        fi
        sleep 10
    done
    wait "$pid"; st=$?
    (( st == 0 )) || { log "export of $i failed ($st, see $x.log)"; rm -f -- "$x.color_sets.txt.zst.part"; return 1; }
    mv -- "$x.color_sets.txt.zst.part" "$x.color_sets.txt.zst" && touch "$x.done"
    log "export of $i done: $(du -h --apparent-size "$x.color_sets.txt.zst" | cut -f1) compressed color sets"
}
# compare_z <description> <prefix A> <prefix B> <index A> <index B>: as check of lib.bash, one line in checks.tsv
compare_z() {
    local d=$1 a=$2 b=$3 ia=$4 ib=$5 c names s v pids=()
    local e=$W/cmp_$$; mkdir -p "$e"
    for s in a b; do
        [[ $s == a ]] && v=$a || v=$b
        ln -sf -- "$v.metadata.txt" "$e/$s.metadata.txt"
        ln -sf -- "$v.unitigs.fa" "$e/$s.unitigs.fa"
        rm -f -- "$e/$s.color_sets.txt"; mkfifo "$e/$s.color_sets.txt"
        nice -n 19 zstd -q -dc "$v.color_sets.txt.zst" > "$e/$s.color_sets.txt" & pids+=($!)
    done
    log "compare: $d"
    if capped "$BIN/compare_color_dumps" "$e/a" "$e/b" -t "$THREADS" > "$e/compare.out" 2>&1; then c=equal; else c=DIFFERENT; fi
    kill "${pids[@]}" 2> /dev/null; wait "${pids[@]}" 2> /dev/null
    if cmp -s <(capped "$BIN/themisto2" dump-color-names -i "$ia" 2> /dev/null) <(capped "$BIN/themisto2" dump-color-names -i "$ib" 2> /dev/null); then names=equal; else names=DIFFERENT; fi
    [[ -s $CHECKS ]] || printf 'check\tkmers_and_color_sets\tcolor_names\tdate\n' > "$CHECKS"
    printf '%s\t%s\t%s\t%s\n' "$d" "$c" "$names" "$(date -Is)" | tee -a "$CHECKS"
    grep -v '^\[' "$e/compare.out" | tail -4
    cp -- "$e/compare.out" "$W/compare_$(echo "$d" | tr -c 'A-Za-z0-9_\n' '_').out"
    rm -f -- "$e"/{a,b}.{metadata.txt,unitigs.fa,color_sets.txt} "$e/compare.out"; rmdir -- "$e"
}
# rm_export <prefix>: the files of export_z, by name
rm_export() { local f; for f in "$1".{metadata.txt,unitigs.fa,color_sets.txt.zst,color_sets.txt.zst.part,log,done}; do [[ -f $f ]] && rm -- "$f"; done; return 0; }
# pair <description> <index A> <index B>
pair() {
    local d=$1 a=$2 b=$3
    done_check "$d" && { log "skip (in checks.tsv): $d"; return 0; }
    [[ -f $a && -f $b ]] || { log "no check of '$d': $a or $b missing"; return 1; }
    ( flock 9
      export_z "$a" "$X/$(basename "$a" .thm2)" && export_z "$b" "$X/$(basename "$b" .thm2)" &&
          compare_z "$d" "$X/$(basename "$a" .thm2)" "$X/$(basename "$b" .thm2)" "$a" "$b"
    ) 9> "$LOCK"
}

if [[ -n ${PAIRS:-} ]]; then
    while IFS='|' read -r d a b; do [[ -n $d ]] && pair "$d" "$a" "$b"; done <<< "$PAIRS"
    log "stream_check done"; exit 0
fi
(( $(free_gb) >= MIN_FREE )) || { log "ABORT: under $MIN_FREE GiB free"; exit 1; }
[[ -f $REBUILD ]] || { log "ABORT: no $REBUILD"; exit 1; }
for n in 121033 60517; do
    out=$W/merge_first_$n.thm2
    if (( MERGES )) && [[ ! -f $out ]]; then
        done_run "pneumo_merge_first_${n}_r$ROUND" && { log "no merge first_$n: round $ROUND is in the ledger but $out is gone (set ROUND)"; continue; }
        [[ -f $PI/first_$n.thm2 && -f $PI/rest_$n.thm2 ]] || { log "no merge first_$n: its inputs are gone"; continue; }
        list=$W/merge_first_$n.list; printf '%s\n' "$PI/first_$n.thm2" "$PI/rest_$n.thm2" > "$list"
        CAP_GB=$CAP INPUTS="$PI/first_$n.thm2 $PI/rest_$n.thm2" run "pneumo_merge_first_${n}_r$ROUND" "$out" 0 \
            "$BIN/themisto2" merge --index-file-list "$list" --temp-dir "$W/tmp" -o "$out" -t "$THREADS" --no-lcs ||
            { log "FAILED merge first_$n"; continue; }
    fi
    d="pneumo merge first_$n = full rebuild (all.direct.thm2), whole index"
    pair "$d" "$out" "$REBUILD" && done_check "$d" && {
        rm_export "$X/merge_first_$n"
        (( ${KEEP_MERGES:-0} )) || { [[ -f $out ]] && rm -- "$out"; }
    }
done
done_check "pneumo merge first_121033 = full rebuild (all.direct.thm2), whole index" &&
    done_check "pneumo merge first_60517 = full rebuild (all.direct.thm2), whole index" && rm_export "$X/all.direct"
log "pneumo subset compl_1 vs its rebuild: not done, the rebuild pneumo_direct_compl_1 was deleted after its n_kmers check (30 Oct)"
log "stream_check done"
