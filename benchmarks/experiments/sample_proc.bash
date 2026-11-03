#!/usr/bin/env bash
# Samples a process every INTERVAL seconds until it exits, from /proc: resident memory (current
# and peak), CPU time of all its threads, and bytes read and written by system calls (rchar and wchar:
# read_bytes stays 0 on network filesystems such as BeeGFS).
# Usage: sample_proc.bash <pid> <out.tsv> [interval_s]
set -u
pid=$1 out=$2 interval=${3:-1}
clk=$(getconf CLK_TCK)
printf 'epoch\trss_kb\thwm_kb\tcpu_s\trchar\twchar\tthreads\n' > "$out"
while [[ -d /proc/$pid ]]; do
    now=$(date +%s.%N)
    status=$(cat "/proc/$pid/status" 2> /dev/null) || break
    stat=$(cat "/proc/$pid/stat" 2> /dev/null) || break
    io=$(cat "/proc/$pid/io" 2> /dev/null) || io=""
    rss=$(awk '/^VmRSS:/ { print $2 }' <<< "$status")
    hwm=$(awk '/^VmHWM:/ { print $2 }' <<< "$status")
    threads=$(awk '/^Threads:/ { print $2 }' <<< "$status")
    # Fields 14 and 15 of /proc/<pid>/stat (utime, stime) count after the command name in parentheses
    cpu=$(awk -v clk="$clk" '{ sub(/^.*\) /, ""); printf "%.2f", ($12 + $13) / clk }' <<< "$stat")
    rd=$(awk '/^rchar:/ { print $2 }' <<< "$io")
    wr=$(awk '/^wchar:/ { print $2 }' <<< "$io")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "${rss:-}" "${hwm:-}" "$cpu" "${rd:-}" "${wr:-}" "${threads:-}" >> "$out"
    sleep "$interval"
done
