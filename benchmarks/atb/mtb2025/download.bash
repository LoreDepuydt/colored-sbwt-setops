#!/bin/bash
# Downloads the ATB assembly archives that hold the Mtb samples of 2024-08 and 2025-05 (lists/archives.tsv,
# from make_lists.py) into archives/, with md5 check. Skips archives already there with the right md5;
# resumes partial downloads. Light on CPU and memory: fine on any server.
# Usage: bash download.bash   (ONLY=<archive name>: just that one)
set -euo pipefail
W=$(cd "$(dirname "$0")" && pwd)
cd "$W/archives"
fail=0
while IFS=$'\t' read -r a url md5 mb _; do
    [[ $a == archive ]] && continue
    [[ -n ${ONLY:-} && $a != "$ONLY" ]] && continue
    if [[ -s $a ]] && [[ $(md5sum < "$a" | cut -d' ' -f1) == "$md5" ]]; then echo "ok   $a (have it)"; continue; fi
    echo "get  $a ($mb MB)"
    for try in 1 2 3; do
        wget -q -c -O "$a.part" "$url" && [[ $(md5sum < "$a.part" | cut -d' ' -f1) == "$md5" ]] && { mv "$a.part" "$a"; break; }
        echo "     retry $try: $a"; [[ $try == 2 ]] && rm -f "$a.part"; sleep 10
    done
    [[ -s $a ]] || { echo "FAILED $a"; fail=1; }
done < "$W/lists/archives.tsv"
n=$(ls "$W"/archives/*.tar.xz 2> /dev/null | wc -l)
echo "$n of $(( $(wc -l < "$W/lists/archives.tsv") - 1 )) archives; $(du -sh "$W/archives" | cut -f1)"
exit $fail
