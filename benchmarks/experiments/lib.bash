# Shared by the scripts of the paper benchmark (source it). Expects P (this folder) to be set, and
# for run/check: THREADS, CAP_GB, W (work folder of the script, under $B/work);
# for run, optionally INPUTS (files read into the page cache before the run, see run.bash).
BIN=${BIN:-$P/bin}
B=${BENCH_DIR:-/tmp/colored-sbwt-bench}
K=31
log() { echo "[$(date -Is)] $*"; }
done_run() { local l=${LEDGER:-$P/ledger.tsv}; [[ -f $l ]] && awk -F'\t' -v n="$1" '$1 == n && $4 == 0 { f = 1 } END { exit !f }' "$l"; }
# run <name> <outputs> <delete outputs 0/1> <command...>: through run.bash, unless done already
run() {
    local name=$1 outputs=$2 del=$3; shift 3
    done_run "$name" && { log "skip $name (in the ledger)"; return 0; }
    # A failed earlier attempt of this run (not in the ledger with exit 0) is set aside, not deleted
    local old=${RESULTS:-$P/runs}/${name}_$(hostname -s)
    [[ -e $old ]] && mv "$old" "$old.failed_$(date +%Y%m%d-%H%M%S)"
    log "run $name"
    mkdir -p "$W/tmp"
    NAME=$name MEM_LIMIT_GB=$CAP_GB THREADS=$THREADS INPUTS=${INPUTS:-} OUTPUTS=$outputs DELETE_OUTPUTS=$del WORKDIR=$W/tmp \
        bash "$P/run.bash" "$@"
}
# guarded <command...>: runs the command, stopped if the disk of $B falls under CHECK_MIN_FREE_GB (300)
# GiB free (an export of a whole S. pneumoniae index wrote 829 GB and filled the disk)
guarded() {
    local pid min=${CHECK_MIN_FREE_GB:-300}
    "$@" & pid=$!
    while kill -0 "$pid" 2> /dev/null; do
        if (( $(df --output=avail -B1G "$B" | tail -1) < min )); then
            kill -TERM "$pid"; wait "$pid"; log "stopped: under $min GiB free on the disk of $B"; return 1
        fi
        sleep 5
    done
    wait "$pid"
}
# check_clean <check folder>: deletes the files check writes, by name (no recursive delete), then the
# folder if it is empty; only for a folder work/<name>/check under $B
check_clean() {
    local e=$1 f
    [[ $e == "$B"/work/?*/check ]] || { log "not deleting $e: not a check folder"; return 1; }
    for f in "$e"/{a,b}.{metadata.txt,unitigs.fa,color_sets.txt,log} "$e/compare.out"; do
        [[ -f $f ]] && rm -- "$f"
    done
    rmdir -- "$e" 2> /dev/null || log "left $e: it holds other files"
}
# check <checks file> <description> <index1> <index2>: exports both (not timed, capped at CAP_GB, stopped
# if the disk runs low) and compares the k-mers and the color set of every k-mer (compare_color_dumps),
# and the color names
check() {
    local checks=$1 d=$2 a=$3 b=$4 e=$W/check c names
    local capped=(systemd-run --user --scope --quiet -p MemoryMax="${CAP_GB}G" -p MemorySwapMax=0 nice -n 19)
    mkdir -p "$e"
    [[ -s $checks ]] || printf 'check\tkmers_and_color_sets\tcolor_names\tdate\n' > "$checks"
    guarded "${capped[@]}" "$BIN/themisto2" export -i "$a" -o "$e/a" -t "$THREADS" 2> "$e/a.log" &&
        guarded "${capped[@]}" "$BIN/themisto2" export -i "$b" -o "$e/b" -t "$THREADS" 2> "$e/b.log" ||
        { log "export failed: $d"; check_clean "$e"; return 1; }
    if "${capped[@]}" "$BIN/compare_color_dumps" "$e/a" "$e/b" -t "$THREADS" > "$e/compare.out" 2>&1; then c=equal; else c=DIFFERENT; fi
    if cmp -s <("${capped[@]}" "$BIN/themisto2" dump-color-names -i "$a" 2> /dev/null) <("${capped[@]}" "$BIN/themisto2" dump-color-names -i "$b" 2> /dev/null); then names=equal; else names=DIFFERENT; fi
    printf '%s\t%s\t%s\t%s\n' "$d" "$c" "$names" "$(date -Is)" | tee -a "$checks"
    tail -3 "$e/compare.out"
    check_clean "$e"
}
