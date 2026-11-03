#!/bin/bash
# Runs one command with the protections of ../experiments/run.bash, without its samples and statistics:
#   - starts only if the server has MEM_LIMIT_GB + MIN_FREE_GB of available memory; otherwise
#     checks again every minute, for at most WAIT_S seconds (default 7200) in all, then gives up
#   - one job at a time: then waits for the lock ../experiments/.lock (shared with the scripts
#     there), held until the command ends; if the memory is gone once it has the lock, it
#     gives the lock back and waits for the memory again (so waiting never holds the lock)
#   - in a cgroup of the user systemd instance capped at MEM_LIMIT_GB, without swap (checked first),
#     at nice 19, a low I/O priority, and oom_score_adj 1000 (the first the kernel kills)
#   - killed by a watchdog if the memory available on the server falls under MIN_FREE_GB, or if the
#     memory pressure goes over MAX_PRESSURE while the kernel has less than PRESSURE_FREE_GB
#     available (another user's processes would otherwise be pushed to swap, which ours cannot use)
# Usage: MEM_LIMIT_GB=<N> bash capped_run.bash <command> [args...]
# Settings: MEM_LIMIT_GB (required), MIN_FREE_GB (32), MAX_PRESSURE (10), PRESSURE_FREE_GB (32),
# WAIT_S (7200). Exit status: that of the command; 125 if it did not start; 137 if the watchdog
# killed it (the reason on stderr, with the 3 largest processes of other users).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
MEM_LIMIT_GB=${MEM_LIMIT_GB:?set MEM_LIMIT_GB}
MIN_FREE_GB=${MIN_FREE_GB:-32}
MAX_PRESSURE=${MAX_PRESSURE:-10}
PRESSURE_FREE_GB=${PRESSURE_FREE_GB:-32}
WAIT_S=${WAIT_S:-7200}
(( $# > 0 )) || { echo "usage: MEM_LIMIT_GB=<N> capped_run.bash <command> [args...]" >&2; exit 125; }
die() { echo "capped_run ABORT: $*" >&2; exit 125; }

# As in ../experiments/run.bash: MemAvailable, plus the part of the ZFS ARC above its minimum (if any)
mem_available_kb() {
    awk '/^MemAvailable:/ { print $2 }' /proc/meminfo |
        awk -v arc=/proc/spl/kstat/zfs/arcstats '{ a = $1 }
            END { while ((getline l < arc) > 0) { split(l, f, " "); if (f[1] == "size") s = f[3]; if (f[1] == "c_min") m = f[3] }
                  print a + (s > m ? int((s - m) / 1024) : 0) }'
}
kernel_available_kb() { awk '/^MemAvailable:/ { print $2 }' /proc/meminfo; }
memory_pressure() { awk '/^full/ { split($2, f, "="); print int(f[2]) }' /proc/pressure/memory; }
others_top() { ps -eo user:32,pid,rss,pcpu,etime,args --sort=-rss --no-headers | awk -v me="$(id -un)" '$1 != me' | head -"$1"; }

need=$(( MEM_LIMIT_GB + MIN_FREE_GB ))
available_gb() { echo $(( $(mem_available_kb) / 1048576 )); }
LOCK=$HERE/../experiments/.lock
exec 9> "$LOCK"
waited=0 told_memory=0 told_lock=0
while :; do
    # The memory first, without the lock, so that a job waiting for memory never holds up the
    # others; the lock has no time limit, the memory WAIT_S in all
    until (( $(available_gb) >= need )); do
        (( waited >= WAIT_S )) && die "only $(available_gb) GiB available on $(hostname) after ${waited} s, need $need"
        (( told_memory++ == 0 )) && echo "capped_run: $(available_gb) GiB available, waiting for $need" >&2
        sleep 60; waited=$(( waited + 60 ))
    done
    flock -n 9 || { (( told_lock++ == 0 )) && echo "capped_run: waiting for the lock $LOCK since $(date -Is)" >&2; flock 9; }
    (( $(available_gb) >= need )) && break
    # The job that held the lock may have left less memory: give the lock back and wait again
    flock -u 9
done

mem_max=$(( MEM_LIMIT_GB << 30 ))
cap=(systemd-run --user --scope --quiet -p MemoryMax="$mem_max" -p MemorySwapMax=0)
probe=$("${cap[@]}" sh -c 'cat "/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)/memory.max"' 2>&1)
[[ $probe == "$mem_max" ]] || die "no memory cap: systemd-run --user --scope gave '$probe', expected $mem_max"
echo "capped_run: on $(hostname) at $(date -Is), cap ${MEM_LIMIT_GB}G, min free ${MIN_FREE_GB}G: $*" >&2

# Every command execs the next one, so $! is the pid of the command itself
"${cap[@]}" choom -n 1000 -- nice -n 19 ionice -c 2 -n 7 "$@" &
pid=$!
trap 'kill -TERM "$pid" 2> /dev/null; exit 130' INT TERM HUP
killed=0
while kill -0 "$pid" 2> /dev/null; do
    reason=""
    if (( $(mem_available_kb) < MIN_FREE_GB * 1048576 )); then
        reason="server available memory under MIN_FREE_GB=$MIN_FREE_GB"
    elif (( $(memory_pressure) > MAX_PRESSURE && $(kernel_available_kb) < PRESSURE_FREE_GB * 1048576 )); then
        reason="memory pressure (full avg10) $(memory_pressure) % over MAX_PRESSURE=$MAX_PRESSURE with $(( $(kernel_available_kb) / 1048576 )) GiB available"
    fi
    if [[ -n $reason ]]; then
        echo "capped_run watchdog: $reason at $(date -Is); killing $*" >&2
        others_top 3 | awk '{ cmd = $6; for (i = 7; i <= NF; i++) cmd = cmd " " $i; printf "  other: %s pid %s, %.1f GiB: %s\n", $1, $2, $3 / 1048576, substr(cmd, 1, 120) }' >&2
        kill -TERM "$pid" 2> /dev/null; sleep 2; kill -KILL "$pid" 2> /dev/null
        killed=1
        break
    fi
    sleep 2
done
wait "$pid"
status=$?
(( killed )) && exit 137
exit $status
