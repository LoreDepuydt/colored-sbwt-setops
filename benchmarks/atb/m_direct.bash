#!/bin/bash
# The subset of the ATB species index to M. tuberculosis (the Mtb colour alone, no --keep-all-colors),
# with the binary of X. As X (mtb_case.bash step 0): capped_run.bash with cap X_CAP (410) GiB and MIN_FREE_GB
# 24, -t 32 --no-lcs, /usr/bin/time -v; without --low-ram-mode first, and only if that run is killed or
# fails, again with --low-ram-mode (logs/M_atb.log, then logs/M_atb_lowram.log; the failed log is kept).
# Result idx/M_atb.thm2, then its stats (idx/M_atb.stats.tsv), which must equal those of M = subset(X, Mtb).
# The mode it ended up using is in logs/M_atb.mode.txt.
# Usage: W=... A=... bash m_direct.bash   (W and A as in mtb_case.bash, after its step 1; it waits up to
#   WAIT_S (2 days) for 434 GiB of available memory)
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
W=${W:?set W, the work folder of mtb_case.bash}
B=${B:-$HERE/../../colors/target/release/themisto2}
CAPPED=$HERE/capped_run.bash
A=${A:?set A, the ATB species index}
MTB=per_species_unitigs/mycobacterium_tuberculosis-unitigs-k31.fna
export WAIT_S=${WAIT_S:-172800} MIN_FREE_GB=${MIN_FREE_GB:-24}
X_CAP=${X_CAP:-410}
log() { echo "[$(date -Is)] $*"; }
cd "$W"
[[ -e idx/M_atb.thm2 ]] && { log "idx/M_atb.thm2 exists: nothing to do"; exit 0; }
md5sum "$B" | tee logs/M_atb.binary.md5
attempt() {  # attempt <log name> [--low-ram-mode]
    local name=$1; shift
    log "run $name"
    MEM_LIMIT_GB=$X_CAP bash "$CAPPED" /usr/bin/time -v "$B" subset --index "$A" -o idx/M_atb.part.thm2 --colors "$MTB" -t 32 --no-lcs "$@" > "logs/$name.log" 2>&1
    local st=$?
    grep -E 'Elapsed|Maximum resident|Exit status' "logs/$name.log" | sed 's/^\s*/    /'
    log "$name: exit $st"
    return $st
}
if attempt M_atb; then
    mode="without --low-ram-mode (logs/M_atb.log)"
else
    rm -f -- idx/M_atb.part.thm2
    attempt M_atb_lowram --low-ram-mode || { rm -f -- idx/M_atb.part.thm2; log "FAILED twice"; exit 1; }
    mode="with --low-ram-mode (logs/M_atb_lowram.log; the attempt without it failed, logs/M_atb.log)"
fi
echo "$mode" > logs/M_atb.mode.txt
log "mode: $mode"
mv -- idx/M_atb.part.thm2 idx/M_atb.thm2
MEM_LIMIT_GB=64 bash "$CAPPED" "$B" stats -i idx/M_atb.thm2 -t 8 > logs/stats_M_atb.log 2>&1
grep -P '^\w+\t' logs/stats_M_atb.log > idx/M_atb.stats.tsv
if diff <(cut -f1,2 idx/M.stats.tsv) <(cut -f1,2 idx/M_atb.stats.tsv) > /dev/null; then log "stats of M_atb = stats of M"; else log "stats DIFFER:"; diff idx/M.stats.tsv idx/M_atb.stats.tsv; fi
log "m_direct done"
