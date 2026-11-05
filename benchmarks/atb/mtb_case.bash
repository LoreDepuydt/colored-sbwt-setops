#!/bin/bash
# Case study, ATB only (plus the human reference CHM13 and two MTBC reference genomes): what in the
# ATB-species colour of M. tuberculosis is the species, and what is shared with other species?
#
#   0+1  X       = subset --keep-all-colors (ATB, Mtb)   Mtb's k-mers, each with all ATB species that have it
#                                                       (the only ATB-scale step)
#   1    M       = subset(X, Mtb)                        the Mtb colour alone; must equal the earlier subset(ATB, Mtb)
#   2    O       = subset(X, all but Mtb)                Mtb k-mers that other species have too, with those species
#        M_priv  = difference(M, O)                      Mtb k-mers in no other ATB species;  |O| + |M_priv| = |M|
#   3    k-mers of O per species, with taxonomy and sylph co-detection in Mtb samples (ATB metadata)
#   4    M_ref   = M ∩ REF,  O_ref = O ∩ REF             where the real Mtb genome is (REF: H37Rv + M. bovis AF2122/97)
#   5    Mp_hum  = M_priv ∩ human,  O_hum = O ∩ human    human sequence, only in Mtb or in many species
#   6    M_rest  = (M_priv − human) − REF                Mtb-only, non-human, non-reference; graph statistics (GFA)
#   9    M_q = subset(X, Mtb) and M_clean = (M_q − O) − human, both queryable: simulated reads of
#        contaminant species and of Mtb (control) against both
#   10   as 2, 4, 5 and 9, but against species only: the colours 'unknown' (73,159 samples without a
#        species call) and the 80 mixed calls ('A;B') are no species, and 'unknown' holds Mtb samples
#        (in step 9, O took the whole Mtb genome with it). O_sp = subset(X, all but Mtb, unknown, mixed);
#        M_priv_sp = M − O_sp; O_sp_ref, O_sp_hum; M_clean2 = (M_q − O_sp) − human, and the reads
#   11   F = the long, low-GC unitigs of the M_rest graph (> F_MIN_KMERS k-mers, GC < F_MAX_GC %): foreign
#        DNA no ATB colour, human or the reference explains. M_clean3 = (M_clean2 − F) ∪ the reference k-mers
#        of M (which O_sp took out with the conserved genes). Checks: |F| + |M_clean2 − F| = |M_clean2|, the
#        reference is complete, human; GC histograms of M_clean2 and M_clean3 (GFA, deleted after gfastat)
#
# Every heavy command goes through capped_run.bash, which waits by itself (up to WAIT_S) for the shared
# lock (taking turns with ../experiments) and for the memory it needs. A step is skipped
# when its output exists: start again after a failure (delete the broken output first).
#
# Inputs (settings): W the work folder; A the ATB species index (contact the authors for a link);
# REF the index of H37Rv (NC_000962) and
# M. bovis AF2122/97 (NC_002945), one color; HUMAN the index of CHM13 v2.0, one color, and CHM13 its
# FASTA (for the human reads); ATB_META the ATB metadata (see mtb_case.py); B the themisto2 binary
# (default: colors/ of this repository); GS gfastat (default: gfastat/ here, cargo build --release);
# OLD_M optional, an earlier subset(ATB, Mtb) that M is compared with.
# Usage: W=... A=... REF=... HUMAN=... CHM13=... ATB_META=... bash mtb_case.bash   (DRY=1: print the
#   commands only; FROM=/TO=<step>: run only those steps)
set -euo pipefail
shopt -s nullglob

HERE=$(cd "$(dirname "$0")" && pwd)
W=${W:?set W, the work folder}
B=${B:-$HERE/../../colors/target/release/themisto2}
GS=${GS:-$HERE/gfastat/target/release/gfastat}
CAPPED=$HERE/capped_run.bash
PY="python3 -I $HERE/mtb_case.py"
A=${A:?set A, the ATB species index (contact the authors for a link)}
OLD_M=${OLD_M:-}
REF=${REF:?set REF, the index of the two reference genomes}
HUMAN=${HUMAN:?set HUMAN, the index of CHM13}
CHM13=${CHM13:?set CHM13, the FASTA of CHM13 v2.0}
export ATB_META=${ATB_META:?set ATB_META, the folder of the ATB metadata}
MTB=${MTB:-per_species_unitigs/mycobacterium_tuberculosis-unitigs-k31.fna}
T=${THREADS:-32}
DRY=${DRY:-0}
FROM=${FROM:-0}
TO=${TO:-11}
NO_LOCK=${NO_LOCK:-0}                              # tests only: no lock, no watchdog
MAX_EXPORT_KMERS=${MAX_EXPORT_KMERS:-6000000000}
F_MIN_KMERS=${F_MIN_KMERS:-301}                    # step 11: unitigs of at least this many k-mers ...
F_MAX_GC=${F_MAX_GC:-60}                           # ... and a GC below this (%) go into F
export WAIT_S=${WAIT_S:-172800}                    # capped_run: wait up to 2 days for the memory
export MIN_FREE_GB=${MIN_FREE_GB:-24}              # capped_run: start with cap + this available, watchdog under it
X_CAP=${X_CAP:-410}                                # step 0; the three earlier ATB subsets peaked at 384-395 GiB
# Step 0 needs --low-ram-mode: without an LCS, the column inversion of the DBG of ATB otherwise takes
# 16 bits per node (~360 GiB on top of the loaded index, over the 410 GiB cap)

cd "$W"
mkdir -p idx exp logs reads res tmp

log() { echo "[$(date -Is)] $*"; }
run() {  # run <cap GiB> <name> <command...>
    local cap=$1 name=$2; shift 2
    log "run $name (cap ${cap}G): $*"
    [[ $DRY == 1 ]] && return 0
    local pre=(env MEM_LIMIT_GB="$cap" bash "$CAPPED")
    [[ $NO_LOCK == 1 ]] && pre=(systemd-run --user --scope -q -p MemoryMax="${cap}G" -p MemorySwapMax=0 nice -n 19)
    if ! "${pre[@]}" /usr/bin/time -v "$@" > "logs/$name.log" 2>&1; then
        log "FAILED $name: see $W/logs/$name.log"; tail -5 "logs/$name.log"; exit 1
    fi
    grep -E 'Elapsed|Maximum resident' "logs/$name.log" | sed 's/^\s*/    /'
}
have() { [[ -s $1 ]] && { log "skip $1 (exists)"; return 0; }; return 1; }
step() { (( FROM <= $1 && $1 <= TO )); }
stats() {
    have "${1%.thm2}.stats.tsv" && return 0
    run 64 "stats_$(basename "${1%.thm2}")" "$B" stats -i "$1" -t 8
    [[ $DRY == 1 ]] || grep -P '^\w+\t' "logs/stats_$(basename "${1%.thm2}").log" > "${1%.thm2}.stats.tsv"
}
nk() { awk '$1=="n_kmers"{print $2}' "${1%.thm2}.stats.tsv"; }
op() {  # op <cap> <intersect|difference> <in1> <in2> <out> [extra]: idx/<out>.thm2, without LCS
    local cap=$1 o=$2 a=$3 b=$4 out=$5; shift 5
    have "idx/$out.thm2" || run "$cap" "$out" "$B" "$o" --index1 "$a" --index2 "$b" -o "idx/$out.thm2" -t "$T" --no-lcs "$@"
    stats "idx/$out.thm2"
}
names() {
    local out=exp/$(basename "${1%.thm2}").names.tsv
    have "$out" && return 0
    log "dump-color-names $1"; [[ $DRY == 1 ]] || "$B" dump-color-names -i "$1" > "$out" 2>/dev/null
}
gfa() {  # gfa <index> [colors]: export as GFA (with colour sets if 'colors') and gfastat
    local name; name=$(basename "${1%.thm2}")
    have "exp/$name.summary.txt" && return 0
    stats "$1"
    if [[ $DRY != 1 ]] && (( $(nk "$1") > MAX_EXPORT_KMERS )); then
        log "REFUSED: $1 has $(nk "$1") k-mers (> MAX_EXPORT_KMERS = $MAX_EXPORT_KMERS); not exporting"; exit 1
    fi
    if [[ ${2:-} == colors ]]; then
        run 96 "export_$name" "$B" export -i "$1" -o "exp/$name" --gfa -t 16
    else
        run 96 "export_$name" "$B" export -i "$1" -o "exp/$name" --gfa --no-color-sets -t 16
    fi
    run 24 "gfastat_$name" "$GS" "exp/$name.unitigs.gfa" "exp/$name"
}
check_sum() {  # check_sum <expected> <index> ...: the n_kmers of the indexes add up to <expected>
    [[ $DRY == 1 ]] && return 0
    local want=$1 total=0 f; shift
    for f in "$@"; do total=$(( total + $(nk "$f") )); done
    if (( total == want )); then log "check: $* add up to $want: ok"; else log "CHECK FAILED: $* add up to $total, not $want"; exit 1; fi
}

if step 0; then
    log "== 0+1: X = subset --keep-all-colors (ATB, Mtb)  [ATB-scale]"
    have idx/X.thm2 || run "$X_CAP" X "$B" subset --index "$A" -o idx/X.thm2 --colors "$MTB" --keep-all-colors -t "$T" --no-lcs --low-ram-mode
    stats idx/X.thm2
fi
if step 1; then
    log "== 1: M = subset(X, Mtb)"
    have idx/M.thm2 || run 96 M "$B" subset --index idx/X.thm2 -o idx/M.thm2 --colors "$MTB" -t "$T" --no-lcs
    stats idx/M.thm2
    if [[ $DRY != 1 ]]; then
        (( $(nk idx/M.thm2) == $(nk idx/X.thm2) )) || { log "CHECK FAILED: |M| != |X|"; exit 1; }
        if [[ -n $OLD_M && -s $OLD_M && ! -s res/M_vs_old.txt ]]; then
            { md5sum idx/M.thm2 "$OLD_M" | awk '{print $1}' | uniq | wc -l | awk '{print ($1 == 1 ? "M is byte-identical to the earlier subset(ATB, Mtb)" : "M differs in bytes from the earlier subset(ATB, Mtb): compare their stats")}'; } | tee res/M_vs_old.txt
        fi
    fi
fi
if step 2; then
    log "== 2: O = subset(X, all but Mtb); M_priv = difference(M, O)"
    have idx/O.thm2 || run 96 O "$B" subset --index idx/X.thm2 -o idx/O.thm2 --colors "$MTB" --complement -t "$T" --no-lcs
    stats idx/O.thm2
    op 96 difference idx/M.thm2 idx/O.thm2 M_priv
    [[ $DRY == 1 ]] || check_sum "$(nk idx/M.thm2)" idx/O.thm2 idx/M_priv.thm2
fi
if step 3; then
    log "== 3: which species share k-mers with Mtb"
    gfa idx/O.thm2 colors; names idx/O.thm2
    have res/O.percolour.tsv || { [[ $DRY == 1 ]] || $PY percolour exp/O exp/O.names.tsv res/O.percolour.tsv | tee res/O.percolour.summary.txt; }
    have res/O.species.tsv || { [[ $DRY == 1 ]] || $PY species res/O.percolour.tsv res/O.species.tsv | tee res/O.species.summary.txt; }
fi
if step 4; then
    log "== 4: the reference genomes"
    op 64 intersect idx/M.thm2 "$REF" M_ref
    op 64 intersect idx/O.thm2 "$REF" O_ref
fi
if step 5; then
    log "== 5: human"
    op 64 intersect  idx/M_priv.thm2 "$HUMAN" Mp_hum
    op 64 intersect  idx/O.thm2      "$HUMAN" O_hum
fi
if step 6; then
    log "== 6: the rest, and the graph statistics"
    op 64 difference idx/M_priv.thm2    "$HUMAN" Mp_nohum
    op 64 difference idx/Mp_nohum.thm2  "$REF"   M_rest
    for p in M_rest Mp_hum M_ref; do gfa idx/$p.thm2; done
fi
if step 9; then
    log "== 9: cleaned index, and simulated reads"
    # queryable versions: without --no-lcs the results get an LCS array
    have idx/M_q.thm2       || run 160 M_q "$B" subset --index idx/X.thm2 -o idx/M_q.thm2 --colors "$MTB" -t "$T"
    have idx/M_q_noO.thm2   || run 96 M_q_noO "$B" difference --index1 idx/M_q.thm2 --index2 idx/O.thm2 -o idx/M_q_noO.thm2 -t "$T"
    have idx/M_clean.thm2   || run 96 M_clean "$B" difference --index1 idx/M_q_noO.thm2 --index2 "$HUMAN" -o idx/M_clean.thm2 -t "$T"
    for x in M_q M_clean; do stats idx/$x.thm2; done
    declare -A REFS=([mtb_H37Rv]=NC_000962 [m_abscessus_ATCC19977]=NC_010397 [m_avium_104]=NC_008595
                     [m_intracellulare_ATCC13950]=NC_016946 [ralstonia_pickettii_12J]=NC_010682
                     [cutibacterium_acnes_KPA171202]=NC_006085 [s_aureus_NCTC8325]=NC_007795 [e_coli_MG1655]=NC_000913)
    for r in "${!REFS[@]}"; do
        have reads/$r.ref.fa && continue
        log "download $r (${REFS[$r]})"; [[ $DRY == 1 ]] && continue
        curl -sf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REFS[$r]}&rettype=fasta&retmode=text" > reads/$r.ref.fa.part
        [[ $(head -c1 reads/$r.ref.fa.part) == ">" ]] || { log "download of $r failed"; exit 1; }
        mv reads/$r.ref.fa.part reads/$r.ref.fa; grep '>' reads/$r.ref.fa | sed 's/^/    /'
    done
    have reads/human_chr21.ref.fa || { log "extract CHM13 chr21"; [[ $DRY == 1 ]] || ( set +o pipefail; zcat "$CHM13" | awk '/^>/{if(p)exit; p=($1==">chr21")} p' > reads/human_chr21.ref.fa ); }
    N=${N_READS:-100000}
    : > reads/queries.txt; : > reads/out_M.txt; : > reads/out_clean.txt
    for f in reads/*.ref.fa; do
        r=$(basename "$f" .ref.fa)
        have reads/$r.reads.fa || { [[ $DRY == 1 ]] || $PY simreads "$f" "$N" 150 1 reads/$r.reads.fa; }
        echo "$W/reads/$r.reads.fa" >> reads/queries.txt
        echo "$W/reads/$r.M.txt" >> reads/out_M.txt
        echo "$W/reads/$r.clean.txt" >> reads/out_clean.txt
    done
    for x in M clean; do
        idx=idx/M_q.thm2; [[ $x == clean ]] && idx=idx/M_clean.thm2
        have reads/done_$x || { run 96 reads_$x "$B" threshold-pseudoalign -i "$idx" --query-list reads/queries.txt --query-output-list reads/out_$x.txt -d 0.7 -n all -t "$T" --themisto1-output-format && { [[ $DRY == 1 ]] || touch reads/done_$x; }; }
    done
    if [[ $DRY != 1 ]]; then
        { printf 'reads\tindex\treads_with_hits\tn_reads\tpercent\n'
          for f in reads/*.reads.fa; do r=$(basename "$f" .reads.fa)
              for x in M clean; do $PY fpsummary "$r	$x" reads/$r.$x.txt "$N"; done; done; } | tee res/reads.tsv
    fi
fi
if step 10; then
    log "== 10: against species only (without unknown and the mixed calls)"
    have idx/O_sp.thm2 || run 96 O_sp "$B" subset --index idx/X.thm2 -o idx/O_sp.thm2 --colors-file keep/not_other_species.txt --complement -t "$T" --no-lcs
    stats idx/O_sp.thm2
    op 96 difference idx/M.thm2 idx/O_sp.thm2 M_priv_sp
    [[ $DRY == 1 ]] || check_sum "$(nk idx/M.thm2)" idx/O_sp.thm2 idx/M_priv_sp.thm2
    op 64 intersect idx/O_sp.thm2 "$REF"   O_sp_ref
    op 64 intersect idx/O_sp.thm2 "$HUMAN" O_sp_hum
    have idx/M_q_noOsp.thm2 || run 96 M_q_noOsp "$B" difference --index1 idx/M_q.thm2 --index2 idx/O_sp.thm2 -o idx/M_q_noOsp.thm2 -t "$T"
    have idx/M_clean2.thm2  || run 96 M_clean2 "$B" difference --index1 idx/M_q_noOsp.thm2 --index2 "$HUMAN" -o idx/M_clean2.thm2 -t "$T"
    stats idx/M_clean2.thm2
    : > reads/out_clean2.txt
    while read -r q; do echo "${q%.reads.fa}.clean2.txt" >> reads/out_clean2.txt; done < reads/queries.txt
    have reads/done_clean2 || { run 96 reads_clean2 "$B" threshold-pseudoalign -i idx/M_clean2.thm2 --query-list reads/queries.txt --query-output-list reads/out_clean2.txt -d 0.7 -n all -t "$T" --themisto1-output-format && { [[ $DRY == 1 ]] || touch reads/done_clean2; }; }
    if [[ $DRY != 1 ]]; then
        N=$(grep -c '>' "$(head -1 reads/queries.txt)")
        { printf 'reads\tindex\treads_with_hits\tn_reads\tpercent\n'
          while read -r q; do r=$(basename "$q" .reads.fa); $PY fpsummary "$r	clean2" "${q%.reads.fa}.clean2.txt" "$N"; done < reads/queries.txt; } | tee res/reads_clean2.tsv
    fi
fi
if step 11; then
    log "== 11: remove the long low-GC Mtb-only unitigs, put the reference back"
    # F: from the GFA of M_rest (step 6), the unitigs with >= F_MIN_KMERS k-mers and GC < F_MAX_GC %, as FASTA
    if ! have exp/F.fa; then
        log "filter exp/M_rest.unitigs.gfa: >= $F_MIN_KMERS k-mers, GC < $F_MAX_GC %"
        [[ $DRY == 1 ]] || { nice -n 19 awk -v k=31 -v min="$F_MIN_KMERS" -v gc="$F_MAX_GC" -F'\t' '
            $1 == "S" { n = length($3); if (n - k + 1 < min) next
                        g = gsub(/[GCgc]/, "&", $3); if (100 * g >= gc * n) next
                        print ">" $2; print $3 }' exp/M_rest.unitigs.gfa > exp/F.fa.part && mv exp/F.fa.part exp/F.fa; }
    fi
    [[ $DRY == 1 ]] || printf '%s\tF\n' "$W/exp/F.fa" > exp/F.fof
    have idx/F.thm2 || run 64 F "$B" build -k 31 --file-colors exp/F.fof -o idx/F.thm2 --temp-dir tmp -t "$T"
    stats idx/F.thm2
    op 96 difference idx/M_clean2.thm2 idx/F.thm2 M_clean2_noF
    [[ $DRY == 1 ]] || check_sum "$(nk idx/M_clean2.thm2)" idx/F.thm2 idx/M_clean2_noF.thm2
    # the reference k-mers of M, in the Mtb colour only (M_ref has the colours Mtb and mtbc)
    have idx/M_ref_mtb.thm2 || run 32 M_ref_mtb "$B" subset --index idx/M_ref.thm2 -o idx/M_ref_mtb.thm2 --colors "$MTB" -t "$T" --no-lcs
    stats idx/M_ref_mtb.thm2
    [[ $DRY == 1 ]] || printf '%s\n' "$W/idx/M_clean2_noF.thm2" "$W/idx/M_ref_mtb.thm2" > tmp/clean3.list
    have idx/M_clean3.thm2 || run 96 M_clean3 "$B" merge --index-file-list tmp/clean3.list --merge-shared-colors --temp-dir tmp -o idx/M_clean3.thm2 -t "$T" --no-lcs
    stats idx/M_clean3.thm2
    op 64 intersect idx/M_clean3.thm2 "$REF"   M_clean3_ref
    op 64 intersect idx/M_clean3.thm2 "$HUMAN" M_clean3_hum
    if [[ $DRY != 1 ]]; then
        (( $(nk idx/M_clean3_ref.thm2) == $(nk idx/M_ref.thm2) )) && log "check: M_clean3 has all $(nk idx/M_ref.thm2) reference k-mers of M: ok" \
            || { log "CHECK FAILED: M_clean3 has $(nk idx/M_clean3_ref.thm2) of the $(nk idx/M_ref.thm2) reference k-mers of M"; exit 1; }
    fi
    for x in M_clean2 M_clean3; do
        gfa idx/$x.thm2
        [[ $DRY == 1 ]] || rm -f exp/$x.unitigs.gfa
    done
    if [[ $DRY != 1 ]]; then
        { printf 'index\tn_kmers\tshare_of_M\tgenomes_of_8.8M\n'
          for x in M M_clean2 F M_clean2_noF M_ref_mtb M_clean3 M_clean3_ref M_clean3_hum; do
              awk -v x=$x -v n="$(nk idx/$x.thm2)" -v m="$(nk idx/M.thm2)" 'BEGIN { printf "%s\t%d\t%.4f\t%.0f\n", x, n, n / m, n / 8.8e6 }'; done; } | tee res/clean3.tsv
    fi
fi
log "== n_kmers per index (both strands)"
for f in idx/*.stats.tsv; do printf '%s\t%s\n' "$(basename "${f%.stats.tsv}")" "$(nk "${f%.stats.tsv}.thm2")"; done
log "done; results in $W/res"
