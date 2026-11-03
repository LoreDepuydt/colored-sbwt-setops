#!/bin/bash
# Checks of the S. pneumoniae merges against the full rebuild, the cheap way first. The two S. pneumoniae merges ("add 1 genome",
# first_121033, and "two halves", first_60517) are run again with --no-lcs, through run.bash, as
# pneumo_merge_first_<n>_check (not a round: the medians of the paper do not include them), and each
# result is compared with the full rebuild work/pneumo_rebuild/all.direct.thm2:
#   1. themisto2 stats (the fields of STATS_FIELDS) and md5. If the md5s match, that is the check.
#   2. If not: the color names, in order (the color IDs are the same in both iff these are equal).
#   3. If the md5s differ and the color order is the same: each of the three indexes is exported with
#      the color sets into sha256sum (no export is stored, except the unitigs and metadata, a few GB),
#      all three at once, and the hashes are compared. Equal hashes = equal exports. The unitigs are
#      hashed independently of their order (see unitig_hash), the color sets as written (same color IDs).
#      If the color order differs, nothing more: that case is for stream_check.bash.
# Before that, the rebuild is copied to DEST (persistent; named all.direct.nolcs.thm2: ggcat -> sbwt -> themisto2
# build --no-lcs, not the original s_pneumoniae.thm2), and the md5 of the copy is checked. If the hashes of
# step 3 differ, nothing more is started: stream_check.bash only after asking.
# One line per check in checks.tsv. Everything that loads an index is under the shared lock (run.bash
# takes it for the merges). The merge results are kept in work/pneumo_fullcheck/ (43 GB and 17 GB).
# Step 3 waits (without holding the lock) for STEP3_NEED_GB (160) GiB of available memory, and is stopped if
# the server falls under MIN_FREE_GB (24) GiB; the sort of the unitigs is capped at 16 GiB.
# Usage: bash merge_check.bash   (settings: DEST (the copy of the rebuild; "" = no copy), THREADS (32),
#   CAP (240 for the merges), CHECK_CAP (100 for each check step), STEP3 (1); for a small test: W_DIR,
#   REBUILD, NS (with W_DIR/merge_first_<n>.thm2 there already, so no merge runs), CHECKS)
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-32}
CAP=${CAP:-240}
CHECK_CAP=${CHECK_CAP:-100}
STEP3=${STEP3:-1}
STEP3_NEED_GB=${STEP3_NEED_GB:-160}   # step 3 starts only with this much memory available (3 exports of about 45 GB)
MIN_FREE_GB=${MIN_FREE_GB:-24}        # step 3 is stopped if the server falls under this
WAIT_S=${WAIT_S:-172800}              # wait at most this long for that memory
CHECKS=${CHECKS:-$P/checks.tsv}
source "$P/lib.bash"
W=${W_DIR:-$B/work/pneumo_fullcheck}; mkdir -p "$W/tmp"
PI=$B/indexes/pneumo/merge
REBUILD=${REBUILD:-$B/work/pneumo_rebuild/all.direct.thm2}
read -r -a NS <<< "${NS:-121033 60517}"
DEST=${DEST-}
LOCK=$P/.lock
STATS_FIELDS="n_kmers n_colors sbwt_length n_sampled_positions n_distinct_color_sets n_sparse_color_sets n_dense_color_sets color_set_storage_size_bytes n_unitigs"
capped() { systemd-run --user --scope --quiet -p MemoryMax="${CHECK_CAP}G" -p MemorySwapMax=0 nice -n 19 "$@"; }
done_check() { [[ -f $CHECKS ]] && awk -F'\t' -v d="$1" '$1 == d { f = 1 } END { exit !f }' "$CHECKS"; }
record() {
    [[ -s $CHECKS ]] || printf 'check\tkmers_and_color_sets\tcolor_names\tdate\n' > "$CHECKS"
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date -Is)" | tee -a "$CHECKS"
}
# The memory available on the server (GiB), as capped_run.bash: MemAvailable plus the ZFS ARC above its minimum
available_gb() {
    awk '/^MemAvailable:/ { print $2 }' /proc/meminfo |
        awk -v arc=/proc/spl/kstat/zfs/arcstats '{ a = $1 }
            END { while ((getline l < arc) > 0) { split(l, f, " "); if (f[1] == "size") s = f[3]; if (f[1] == "c_min") m = f[3] }
                  print int((a + (s > m ? int((s - m) / 1024) : 0)) / 1048576) }'
}
field() { awk -F'\t' -v k="$1" '$1 == k { print $2 }' "$2"; }
md5_of() { [[ -s $1.md5 && $1.md5 -nt $1 ]] || nice -n 19 ionice -c3 md5sum "$1" | cut -d' ' -f1 > "$1.md5"; cat "$1.md5"; }
# hash_export <index> <prefix>: sha256 of the color sets (streamed), of the unitigs and of the metadata
hash_export() {
    local h
    for h in colors unitigs metadata; do rm -f -- "$2.$h.sha256" "$2.$h.sha256.part"; done
    capped "$BIN/themisto2" export -i "$1" -o "$2" -t "$THREADS" --colors-to-stdout 2> "$2.log" | sha256sum | cut -d' ' -f1 > "$2.colors.sha256.part" &&
        unitig_hash "$2.unitigs.fa" > "$2.unitigs.sha256.part" && sha256sum < "$2.metadata.txt" | cut -d' ' -f1 > "$2.metadata.sha256.part" ||
        { for h in colors unitigs metadata; do rm -f -- "$2.$h.sha256.part"; done; return 1; }
    # only a complete run leaves hashes (a stopped export would hash a truncated stream)
    for h in colors unitigs metadata; do mv -- "$2.$h.sha256.part" "$2.$h.sha256"; done
}
# unitig_hash <unitigs.fa>: sha256 of the sorted lines "<unitig, in its smaller orientation> <color set id>"
# (export writes the unitigs in a different order on each run; the unitig ids are left out for that reason)
unitig_hash() {
    awk '/^>/ { sub(/.*color_set_id=/, ""); c = $1; next } { print $0 "\t" c }' "$1" > "$1.lines" &&
        paste "$1.lines" <(cut -f1 "$1.lines" | rev | tr ACGT TGCA) |
        awk -F'\t' '{ print ($1 <= $3 ? $1 : $3) "\t" $2 }' |
        systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 nice -n 19 env LC_ALL=C sort -S 8G --parallel=8 -T "$W/tmp" |
        sha256sum | cut -d' ' -f1
    local st=$?; rm -f -- "$1.lines"; return $st
}

[[ -f $REBUILD ]] || { log "ABORT: no $REBUILD"; exit 1; }
# 0. The persistent copy of the rebuild
if [[ -n $DEST && ! -f $DEST ]]; then
    log "copy $REBUILD -> $DEST"
    nice -n 19 ionice -c3 cp -- "$REBUILD" "$DEST.part" && mv -- "$DEST.part" "$DEST" || { log "copy failed"; exit 1; }
fi

# The merges, kept
for n in "${NS[@]}"; do
    out=$W/merge_first_$n.thm2
    [[ -f $out ]] && continue
    done_run "pneumo_merge_first_${n}_check" && { log "no merge first_$n: it is in the ledger but $out is gone"; continue; }
    list=$W/merge_first_$n.list; printf '%s\n' "$PI/first_$n.thm2" "$PI/rest_$n.thm2" > "$list"
    CAP_GB=$CAP INPUTS="$PI/first_$n.thm2 $PI/rest_$n.thm2" run "pneumo_merge_first_${n}_check" "$out" 0 \
        "$BIN/themisto2" merge --index-file-list "$list" --temp-dir "$W/tmp" -o "$out" -t "$THREADS" --no-lcs ||
        log "FAILED merge first_$n"
done

( flock 9
  if [[ -n $DEST ]]; then
      a=$(md5_of "$REBUILD"); b=$(nice -n 19 ionice -c3 md5sum "$DEST" | cut -d' ' -f1)
      log "md5 rebuild $a, copy $b: $( [[ $a == "$b" ]] && echo equal || echo DIFFERENT)"
  fi
  [[ -s $W/rebuild.stats.tsv ]] || capped "$BIN/themisto2" stats -i "$REBUILD" -t "$THREADS" 2> /dev/null > "$W/rebuild.stats.tsv"
  capped "$BIN/themisto2" dump-color-names -i "$REBUILD" 2> /dev/null > "$W/rebuild.color_names.txt"
  : > "$W/step3.txt"
  for n in "${NS[@]}"; do
      out=$W/merge_first_$n.thm2 d="pneumo merge first_$n = full rebuild (all.direct.thm2), whole index"
      done_check "$d" && continue
      [[ -f $out ]] || { log "no check of first_$n: $out missing"; continue; }
      capped "$BIN/themisto2" stats -i "$out" -t "$THREADS" 2> /dev/null > "$W/merge_first_$n.stats.tsv"
      diff=""
      for k in $STATS_FIELDS; do
          x=$(field "$k" "$W/merge_first_$n.stats.tsv"); y=$(field "$k" "$W/rebuild.stats.tsv")
          [[ $x == "$y" ]] || diff+="${diff:+; }$k $x vs $y"
      done
      record "pneumo merge first_$n: themisto2 stats = of the full rebuild" "$( [[ -z $diff ]] && echo equal || echo "DIFFERENT ($diff)")" n/a
      capped "$BIN/themisto2" dump-color-names -i "$out" 2> /dev/null > "$W/merge_first_$n.color_names.txt"
      cmp -s "$W/merge_first_$n.color_names.txt" "$W/rebuild.color_names.txt" && names=equal || names=DIFFERENT
      ma=$(md5_of "$out"); mb=$(md5_of "$REBUILD")
      log "first_$n: md5 $ma, rebuild $mb; color names in order: $names"
      if [[ $ma == "$mb" ]]; then
          record "$d (md5)" equal "$names"
      elif [[ $names == equal ]]; then
          echo "$n" >> "$W/step3.txt"
      else
          log "first_$n: the color order differs: no hash check here, use stream_check.bash"
      fi
  done
  log "md5 and stats done"
) 9> "$LOCK"

# 3. Waits (without the lock) for STEP3_NEED_GB of memory, then, under the lock, the three exports at once,
# each its own process group, all stopped if the server falls under MIN_FREE_GB available
mapfile -t step3 < "$W/step3.txt"
if (( STEP3 && ${#step3[@]} )); then
    waited=0
    while :; do
        until (( $(available_gb) >= STEP3_NEED_GB )); do
            (( waited >= WAIT_S )) && { log "step 3 not started: only $(available_gb) GiB available after $waited s"; exit 1; }
            (( waited % 3600 == 0 )) && log "step 3: $(available_gb) GiB available, waiting for $STEP3_NEED_GB"
            sleep 60; (( waited += 60 ))
        done
        ( flock 9
          (( $(available_gb) >= STEP3_NEED_GB )) || exit 75   # the memory went while waiting for the lock
          set -m
          log "step 3: exports into sha256sum, at once: the rebuild and first_${step3[*]} ($(available_gb) GiB available)"
          pids=()
          [[ -s $W/rebuild.colors.sha256 && -s $W/rebuild.unitigs.sha256 && -s $W/rebuild.metadata.sha256 ]] || { hash_export "$REBUILD" "$W/rebuild" & pids+=($!); }
          for n in "${step3[@]}"; do hash_export "$W/merge_first_$n.thm2" "$W/merge_first_$n" & pids+=($!); done
          stopped=0
          while :; do
              alive=0; for p in "${pids[@]}"; do kill -0 "$p" 2> /dev/null && alive=1; done
              (( alive )) || break
              if (( $(available_gb) < MIN_FREE_GB )); then
                  log "step 3 stopped: under $MIN_FREE_GB GiB available on the server"; stopped=1
                  for p in "${pids[@]}"; do kill -TERM -- "-$p" 2> /dev/null; done
              fi
              sleep 5
          done
          st=0; for p in "${pids[@]}"; do wait "$p" || st=1; done
          (( st == 0 )) || log "step 3: an export failed or was stopped (see $W/*.log)"
          if (( ! stopped )); then
              for n in "${step3[@]}"; do
                  d="pneumo merge first_$n = full rebuild (all.direct.thm2), whole index"
                  [[ -s $W/merge_first_$n.metadata.sha256 && -s $W/rebuild.metadata.sha256 ]] || { log "first_$n: no complete hashes"; continue; }
                  r=equal
                  for h in colors unitigs metadata; do cmp -s "$W/merge_first_$n.$h.sha256" "$W/rebuild.$h.sha256" || r=DIFFERENT; done
                  log "first_$n sha256: $(paste -d' ' "$W"/merge_first_$n.{colors,unitigs,metadata}.sha256) vs rebuild $(paste -d' ' "$W"/rebuild.{colors,unitigs,metadata}.sha256)"
                  record "$d (sha256 of the export)" "$r" equal
                  [[ $r == equal ]] || log "first_$n: the hashes differ; stream_check.bash is NOT started (ask first)"
              done
          fi
          for f in "$W"/rebuild.unitigs.fa "$W"/merge_first_*.unitigs.fa "$W"/*.unitigs.fa.lines "$W"/*.sha256.part; do [[ -f $f ]] && rm -- "$f"; done
          exit 0
        ) 9> "$LOCK"
        (( $? == 75 )) || break
        log "step 3: memory gone once the lock was taken; waiting again"
    done
fi
log "merge_check done"
