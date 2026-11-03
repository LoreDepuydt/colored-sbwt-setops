#!/bin/bash
# Downloads the ATB assemblies (<accession>.fa.gz, kept compressed) of a set of sets/pneumo/ from the
# public bucket allthebacteria-assemblies (eu-west-2, no account needed) into
# ${BENCH_DIR:-/tmp/colored-sbwt-bench}/data/pneumo/genomes/, and checks each against its md5 in sets/pneumo/atb.tsv.
# Files already there with the right md5 are skipped; a file with a wrong md5 is downloaded again,
# and the accessions that still fail are listed at the end. Light: JOBS parallel downloads (default 8),
# at nice 19. Usage: bash download_pneumo.bash <set, e.g. keep_10>
set -uo pipefail
P=$(cd "$(dirname "$0")" && pwd)
SET=${1:?usage: download_pneumo.bash <set>}
JOBS=${JOBS:-8}
OUT=${BENCH_DIR:-/tmp/colored-sbwt-bench}/data/pneumo/genomes
URL=https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com
mkdir -p "$OUT"
[[ -f $P/sets/pneumo/$SET.txt ]] || { echo "no set $P/sets/pneumo/$SET.txt" >&2; exit 1; }
avail=$(df --output=avail -B1G "$OUT" | tail -1 | tr -d ' ')
need=$(awk -F'\t' 'NR == FNR { s[$1] = 1; next } FNR > 1 && ($1 in s) { b += $3 } END { print int(b / 2^30) + 1 }' \
    "$P/sets/pneumo/$SET.txt" "$P/sets/pneumo/atb.tsv")
(( avail >= need + 500 )) || { echo "only $avail GiB free, need $need + 500 GiB margin" >&2; exit 1; }
echo "[$(date -Is)] $SET: $(wc -l < "$P/sets/pneumo/$SET.txt") genomes, ~$need GiB, into $OUT"

# One line per genome of the set: accession md5
awk -F'\t' 'NR == FNR { s[$1] = 1; next } FNR > 1 && ($1 in s) { print $2, $4 }' "$P/sets/pneumo/$SET.txt" "$P/sets/pneumo/atb.tsv" |
    nice -n 19 xargs -P "$JOBS" -n 2 sh -c '
        a=$1 md5=$2 f=$0/$1.fa.gz
        [ -f "$f" ] && [ "$(md5sum < "$f" | cut -d" " -f1)" = "$md5" ] && exit 0
        for try in 1 2 3; do
            curl -sf --retry 3 -o "$f.part" "'"$URL"'/$a.fa.gz" &&
                [ "$(md5sum < "$f.part" | cut -d" " -f1)" = "$md5" ] && mv "$f.part" "$f" && exit 0
            sleep 5
        done
        rm -f "$f.part"; echo "FAILED $a" >&2' "$OUT" 2> "$OUT/../download_$SET.err"
n_ok=$(awk '{ print $1 ".fa.gz" }' "$P/sets/pneumo/$SET.txt" | sed 's#.*/##; s/\.fa\.fa\.gz$/.fa.gz/' | (cd "$OUT" && xargs ls 2> /dev/null) | wc -l)
echo "[$(date -Is)] $SET: $n_ok of $(wc -l < "$P/sets/pneumo/$SET.txt") present; failures: $(grep -c FAILED "$OUT/../download_$SET.err")"
