#!/bin/bash
# The M. tuberculosis colour of ATB 2024-08 and of 2025-05, built from the ATB assemblies, by a merge tree:
#
#   build    per archive (51, lists/archives.tsv): unpack its Mtb assemblies, one GGCAT unitig file per
#            group (K = Mtb in both releases, D = Mtb in 2024-08 only, N = Mtb in 2025-05 only; lists/ from
#            make_lists.py), themisto2 build --from-unitigs with the colours K, D, N; delete the FASTA
#   merge    one multi-way merge of the 51 indexes, --merge-shared-colors: KDN.thm2 (colours K, D, N)
#   compare  Mtb_2024 = subset(KDN, K D), Mtb_2025 = subset(KDN, K N);
#            against the ATB-species colour M (../mtb_atb): M − Mtb_2024 and Mtb_2024 − M (both 0 if the two
#            routes agree); lost = Mtb_2024 − Mtb_2025, new = Mtb_2025 − Mtb_2024; lost ∩ F, O_sp, human and
#            M_clean3 (what the index-only cleaning of ../mtb_atb removed or kept); F and O_sp ∩ Mtb_2025
#   test     one archive (ARCHIVE, default the smallest) through 'build', plus the check of the route: the
#            first TEST_N assemblies of it built from unitigs and directly from the FASTA (plain themisto2
#            build) must give the same k-mers and colour sets (compare_color_dumps)
#
# Every command runs through ../../experiments/run.bash (timed, sampled, cold start, its
# lock shared with the paper benchmark: runs take turns with it and wait for a quiet server), with its
# own ledger (ledger.tsv) and run folders (runs/) here. A run in the ledger with exit 0 is skipped: start
# again after a failure. Heavy files go to the local scratch S; the
# results are kept here.
#
# Usage: MODE=test|build|merge|compare|all bash mtb2025.bash   (settings below)
set -euo pipefail
W=$(cd "$(dirname "$0")" && pwd)
P=$(cd "$W/../../experiments" && pwd)
BIN=${BIN:-$P/bin}                                # the binaries, as for ../../experiments
S=${S:?set S, the local scratch and results}
MA=${MA:?set MA, the idx/ folder of ../mtb_case.bash}   # the 2024-08 case study
HUMAN=${HUMAN:?set HUMAN, the index of CHM13}
MODE=${MODE:?set MODE=test, build, merge, compare or all}
T=${THREADS:-32}
GGCAT_MEM=${GGCAT_MEM:-8}
MERGE_CAP=${MERGE_CAP:-240}
TEST_N=${TEST_N:-50}
export LEDGER=$W/ledger.tsv RESULTS=$W/runs WORK_ROOT=$S
mkdir -p "$S/x" "$S/idx" "$S/tmp" "$S/lists" "$W/res"

log() { echo "[$(date -Is)] $*"; }
done_run() { [[ -f $LEDGER ]] && awk -F'\t' -v n="$1" '$1 == n && $4 == 0 { f = 1 } END { exit !f }' "$LEDGER"; }
run() {  # run <name> <cap GiB> <outputs> <command...>: through run.bash, unless in the ledger
    local name=$1 cap=$2 outputs=$3; shift 3
    done_run "$name" && { log "skip $name (in the ledger)"; return 0; }
    local old=$RESULTS/${name}_$(hostname -s)
    [[ -e $old ]] && mv "$old" "$old.failed_$(date +%Y%m%d-%H%M%S)"
    log "run $name"
    NAME=$name MEM_LIMIT_GB=$cap THREADS=$T INPUTS=${INPUTS:-} OUTPUTS=$outputs DELETE_OUTPUTS=0 WORKDIR=$S/tmp \
        bash "$P/run.bash" "$@" || { log "FAILED $name: see $RESULTS/${name}_$(hostname -s)"; return 1; }
}
capped() { systemd-run --user --scope --quiet -p MemoryMax="${CAP:-64}G" -p MemorySwapMax=0 nice -n 19 "$@"; }
nk() {  # nk <index>: its n_kmers (themisto2 stats, kept next to it)
    local st=${1%.thm2}.stats.tsv
    [[ -s $st && $st -nt $1 ]] || capped "$BIN/themisto2" stats -i "$1" -t 8 2> /dev/null | grep -P '^\w+\t' > "$st"
    awk '$1 == "n_kmers" { print $2 }' "$st"
}
in_scratch_rm() {  # rm -rf, only strictly inside $S
    local p; p=$(realpath -m "$1")
    [[ $p == "$(realpath -m "$S")"/?* ]] || { log "refusing to delete $1: not inside $S"; return 1; }
    rm -rf -- "$p"
}
id_of() { local a=${1#atb.assembly.}; echo "${a%.tar.xz}"; }

# unitig_build <archive> <id> <sample filter: all, or the first N> <output>: unpack, GGCAT per group,
# themisto2 build --from-unitigs (the runs are named <step>_<id>)
unitig_build() {
    local a=$1 id=$2 n=$3 out=$4 x=$S/x/$2 fof=$S/lists/$2.fof g
    done_run "build_$id" && [[ -s $out ]] && { log "skip $id (built)"; return 0; }
    [[ -s $W/archives/$a ]] || { log "missing archive $a: run download.bash"; return 1; }
    awk -F'\t' -v a="$a" -v n="$n" '$3 == a && (n == "all" || ++c <= n)' "$W/lists/samples.tsv" > "$S/lists/$id.samples.tsv"
    cut -f4 "$S/lists/$id.samples.tsv" > "$S/lists/$id.paths"
    in_scratch_rm "$x" 2> /dev/null || true; mkdir -p "$x"
    INPUTS="$W/archives/$a" run "extract_$id" 8 "" tar -xJf "$W/archives/$a" -C "$x" -T "$S/lists/$id.paths" || return 1
    : > "$fof"
    for g in K D N; do
        awk -F'\t' -v g=$g -v x="$x" '$2 == g { print x "/" $4 }' "$S/lists/$id.samples.tsv" > "$S/lists/$id.$g.txt"
        [[ -s $S/lists/$id.$g.txt ]] || continue
        INPUTS="@$S/lists/$id.$g.txt" run "ggcat_${id}_$g" $(( GGCAT_MEM * 3 )) "$x/$g.unitigs.fa" \
            "$BIN/ggcat" build -l "$S/lists/$id.$g.txt" -s 1 -k 31 -m "$GGCAT_MEM" -p -j "$T" -t "$S/tmp" -o "$x/$g.unitigs.fa" || return 1
        printf '%s\t%s\n' "$x/$g.unitigs.fa" "$g" >> "$fof"
    done
    INPUTS="@$fof" run "build_$id" 64 "$out" \
        "$BIN/themisto2" build --from-unitigs --file-colors "$fof" -o "$out" -k 31 -d 30 -t "$T" --temp-dir "$S/tmp" --no-lcs || return 1
    log "$id: $(nk "$out") k-mers"
}

if [[ $MODE == test ]]; then
    a=${ARCHIVE:-$(awk -F'\t' 'NR > 1 { print $4 "\t" $1 }' "$W/lists/archives.tsv" | sort -g | head -1 | cut -f2)}
    id=$(id_of "$a")
    log "== test on $a"
    # 1. the route on the first TEST_N assemblies, against a plain build from the FASTA
    tid=test${TEST_N}_$id
    unitig_build "$a" "$tid" "$TEST_N" "$S/idx/$tid.unitigs.thm2"
    # in the order K, D, N, as the unitig files: the same colour ids
    for g in K D N; do awk -F'\t' -v g=$g -v x="$S/x/$tid" '$2 == g { print x "/" $4 "\t" $2 }' "$S/lists/$tid.samples.tsv"; done > "$S/lists/$tid.plain.fof"
    INPUTS="@$S/lists/$tid.plain.fof" run "plainbuild_$tid" 64 "$S/idx/$tid.plain.thm2" \
        "$BIN/themisto2" build --file-colors "$S/lists/$tid.plain.fof" -o "$S/idx/$tid.plain.thm2" -k 31 -d 30 -t "$T" --temp-dir "$S/tmp" --no-lcs
    e=$S/x/check; mkdir -p "$e"
    CAP=64 capped "$BIN/themisto2" export -i "$S/idx/$tid.unitigs.thm2" -o "$e/a" -t "$T" 2> "$e/a.log"
    CAP=64 capped "$BIN/themisto2" export -i "$S/idx/$tid.plain.thm2" -o "$e/b" -t "$T" 2> "$e/b.log"
    if CAP=64 capped "$BIN/compare_color_dumps" "$e/a" "$e/b" -t "$T" > "$e/compare.out" 2>&1; then c=equal; else c=DIFFERENT; fi
    printf 'route check (%s assemblies of %s): unitigs vs plain build: %s; k-mers %s vs %s\n' "$TEST_N" "$id" "$c" \
        "$(nk "$S/idx/$tid.unitigs.thm2")" "$(nk "$S/idx/$tid.plain.thm2")" | tee "$W/res/test_route_check.txt"
    tail -3 "$e/compare.out"
    in_scratch_rm "$e"; in_scratch_rm "$S/x/$tid"
    # 2. the whole archive, as in 'build' (its index is kept and reused by 'build')
    unitig_build "$a" "$id" all "$S/idx/$id.thm2"
    in_scratch_rm "$S/x/$id"
    awk -F'\t' -v id="$id" 'NR == 1 || $1 ~ id' "$LEDGER" | cut -f1,4,5,8 | column -t
    log "test done"
fi

if [[ $MODE == build || $MODE == all ]]; then
    log "== build: one index per archive"
    while IFS=$'\t' read -r a _; do
        [[ $a == archive ]] && continue
        id=$(id_of "$a")
        unitig_build "$a" "$id" all "$S/idx/$id.thm2" || exit 1
        in_scratch_rm "$S/x/$id"
    done < "$W/lists/archives.tsv"
fi

if [[ $MODE == merge || $MODE == all ]]; then
    log "== merge of the 51 indexes"
    awk -F'\t' -v s="$S/idx" 'NR > 1 { a = $1; sub(/^atb\.assembly\./, "", a); sub(/\.tar\.xz$/, "", a); print s "/" a ".thm2" }' \
        "$W/lists/archives.tsv" > "$S/lists/merge.list"
    while read -r f; do [[ -s $f ]] || { log "missing $f: run MODE=build"; exit 1; }; done < "$S/lists/merge.list"
    INPUTS="@$S/lists/merge.list" run merge_KDN "$MERGE_CAP" "$S/idx/KDN.thm2" \
        "$BIN/themisto2" merge --index-file-list "$S/lists/merge.list" --merge-shared-colors --temp-dir "$S/tmp" -o "$S/idx/KDN.thm2" -t "$T" --no-lcs || exit 1
    log "KDN: $(nk "$S/idx/KDN.thm2") k-mers; colours: $("$BIN/themisto2" dump-color-names -i "$S/idx/KDN.thm2" 2> /dev/null | cut -f2 | tr '\n' ' ')"
fi

if [[ $MODE == compare || $MODE == all ]]; then
    log "== compare with the ATB-species colour and the cleaning"
    I=$S/idx
    op() {  # op <name> <subset|intersect|difference> <args...>: $I/<name>.thm2
        local name=$1 o=$2; shift 2
        if [[ $o == subset ]]; then
            INPUTS="$1" run "$name" 96 "$I/$name.thm2" "$BIN/themisto2" subset --index "$1" --colors "$2" -o "$I/$name.thm2" -t "$T" --no-lcs || exit 1
        else
            INPUTS="$1 $2" run "$name" 96 "$I/$name.thm2" "$BIN/themisto2" "$o" --index1 "$1" --index2 "$2" -o "$I/$name.thm2" -t "$T" --no-lcs || exit 1
        fi
    }
    op Mtb_2024 subset "$I/KDN.thm2" K,D
    op Mtb_2025 subset "$I/KDN.thm2" K,N
    op M_not_built difference "$MA/M.thm2" "$I/Mtb_2024.thm2"
    op built_not_M difference "$I/Mtb_2024.thm2" "$MA/M.thm2"
    op lost difference "$I/Mtb_2024.thm2" "$I/Mtb_2025.thm2"
    op new  difference "$I/Mtb_2025.thm2" "$I/Mtb_2024.thm2"
    op lost_F    intersect "$I/lost.thm2" "$MA/F.thm2"
    op lost_O_sp intersect "$I/lost.thm2" "$MA/O_sp.thm2"
    op lost_hum  intersect "$I/lost.thm2" "$HUMAN"
    [[ -s $MA/M_clean3.thm2 ]] && op lost_kept intersect "$I/lost.thm2" "$MA/M_clean3.thm2"
    op F_in_2025    intersect "$MA/F.thm2" "$I/Mtb_2025.thm2"
    op O_sp_in_2025 intersect "$MA/O_sp.thm2" "$I/Mtb_2025.thm2"
    { printf 'index\tn_kmers\tmeaning\n'
      while IFS=$'\t' read -r f what; do
          [[ -s $f ]] && printf '%s\t%s\t%s\n' "$(basename "${f%.thm2}")" "$(nk "$f")" "$what"
      done <<EOF
$MA/M.thm2	the ATB-species colour of Mtb (2024-08)
$I/Mtb_2024.thm2	built from the assemblies: Mtb in 2024-08
$I/Mtb_2025.thm2	built from the assemblies: Mtb in 2025-05
$I/M_not_built.thm2	in M, not in Mtb_2024 (0 if the routes agree)
$I/built_not_M.thm2	in Mtb_2024, not in M (0 if the routes agree)
$I/lost.thm2	Mtb_2024 − Mtb_2025: only the dropped samples had them
$I/new.thm2	Mtb_2025 − Mtb_2024: only the new samples have them
$I/lost_F.thm2	lost, and in F (removed by the cleaning: long low-GC pieces)
$I/lost_O_sp.thm2	lost, and shared with another ATB species (removed by the cleaning)
$I/lost_hum.thm2	lost, and human
$I/lost_kept.thm2	lost, but kept by the cleaning (in M_clean3)
$I/F_in_2025.thm2	removed by the cleaning (F), but still in Mtb 2025
$I/O_sp_in_2025.thm2	removed by the cleaning (O_sp), but still in Mtb 2025
EOF
    } | tee "$W/res/compare.tsv"
fi
log "done ($MODE)"
