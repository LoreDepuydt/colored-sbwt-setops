#!/bin/bash
# One timed and sampled run of any command (themisto2 set operations, builds, ggcat, sbwt) for the
# paper benchmark, on a shared server.
#
# Before starting, the run waits (checking every minute, for at most WAIT_S seconds, 0 = forever)
#   - for the lock .lock in this folder, shared by all scripts here
#   - for a quiet server: at least MEM_LIMIT_GB + MIN_FREE_GB of available memory, and at least
#     MIN_IDLE_CPUS idle logical CPUs (default THREADS + 8, at most all of them) on average over 30 seconds
# and refuses to start (exit 125) if the local disk has less than MIN_DISK_FREE_GB free, or the work
# folder WORK_ROOT already holds more than DISK_BUDGET_GB.
# Available memory is the kernel's MemAvailable plus the part of the ZFS ARC above its minimum c_min:
# the kernel counts the ARC as used, though it shrinks under memory pressure.
#
# Then the input files INPUTS are dropped from the page cache (written back with sync, then
# dd iflag=nocache, which needs no root), so that every run starts cold and reads its inputs from the
# disk; what is still cached after that (cached_bytes.py, e.g. a file another process keeps open)
# goes to cache.txt.
#
# The command runs in its own cgroup (systemd-run --user --scope, capped at MEM_LIMIT_GB, no swap,
# checked from inside), at nice 19, a low I/O priority and oom_score_adj 1000, under /usr/bin/time -v,
# and is killed by a watchdog if the server runs short of memory: available memory under MIN_FREE_GB;
# the kernel's own available memory (without the ARC) under KERNEL_MIN_FREE_GB for 10 s, when the
# ARC does not shrink fast enough; or all tasks stalled on memory more than MAX_PRESSURE % of the
# last 10 s (/proc/pressure/memory, full avg10) while the kernel has less than PRESSURE_FREE_GB available.
# The watchdog also kills it if the disk of WORK_ROOT falls under RUN_MIN_DISK_FREE_GB free.
#
# Usage: NAME=<run name> MEM_LIMIT_GB=<N> THREADS=<n> [OUTPUTS="<files>"] bash run.bash <command> [args...]
# Settings: MIN_FREE_GB (64), MIN_IDLE_CPUS (THREADS + 8), WAIT_S (0 = wait forever),
#   MIN_DISK_FREE_GB (200), DISK_BUDGET_GB (600), SAMPLE_S (0.2), MAX_PRESSURE (10), PRESSURE_FREE_GB (32),
#   KERNEL_MIN_FREE_GB (4), RUN_MIN_DISK_FREE_GB (100),
#   INPUTS (files or folders the command reads, dropped from the cache before the run; @<file> stands for the paths
#   in the first column of that file, e.g. a list of genomes),
#   OUTPUTS (files the command writes: size and md5 recorded), DELETE_OUTPUTS=1 (then delete them;
#   only files under WORK_ROOT), WORKDIR (a folder under WORK_ROOT the command uses for its temporary
#   files: its size is sampled every 2 s, and it is emptied after the run), STDOUT (where stdout
#   goes, default <run>/stdout.txt), LEDGER (default paper/ledger.tsv), RESULTS (paper/runs), NOTE.
#
# Everything goes to $RESULTS/<NAME>_<host>/ (on NFS, small files only):
#   cmd.txt       the command, settings, host, start time, server state at the start
#   cache.txt     the inputs dropped from the page cache before the run: files, bytes, bytes still cached
#   stderr.log    stderr of the command, each line prefixed with the epoch in microseconds (as ts '%.s')
#   stdout.txt    stdout of the command (unless STDOUT is given)
#   time.txt      /usr/bin/time -v of the command
#   samples.tsv   every SAMPLE_S s: RSS, peak RSS, CPU seconds, bytes read and written, threads
#   server.tsv    every 2 s: memory of the cgroup (all, anonymous, file cache), its CPU (usec), the
#                 CPU of the whole server (busy and total jiffies), available memory (with the ARC),
#                 pressure, the kernel's available memory, the ARC size
#   workdir.tsv   every 2 s: size of WORKDIR in KiB (temporary disk)
#   others.tsv    every 30 s: the 5 processes of other users with the most memory
#   result.txt    one line: the fields of the ledger
# and one line in LEDGER (exit 0 or not). Exit status: that of the command; 125 if it did not start;
# 137 if the watchdog killed it (reason in result.txt and killed.txt).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LOCK=$HERE/.lock
WORK_ROOT=${WORK_ROOT:-${BENCH_DIR:-/tmp/colored-sbwt-bench}/work}
NAME=${NAME:?set NAME}
MEM_LIMIT_GB=${MEM_LIMIT_GB:?set MEM_LIMIT_GB}
THREADS=${THREADS:?set THREADS (the threads the command uses, for the idle CPU check)}
MIN_FREE_GB=${MIN_FREE_GB:-64}
MIN_IDLE_CPUS=${MIN_IDLE_CPUS:-$(( THREADS + 8 < $(nproc) ? THREADS + 8 : $(nproc) ))}
WAIT_S=${WAIT_S:-0}
MIN_DISK_FREE_GB=${MIN_DISK_FREE_GB:-200}
DISK_BUDGET_GB=${DISK_BUDGET_GB:-600}
SAMPLE_S=${SAMPLE_S:-0.2}
MAX_PRESSURE=${MAX_PRESSURE:-10}
PRESSURE_FREE_GB=${PRESSURE_FREE_GB:-32}
KERNEL_MIN_FREE_GB=${KERNEL_MIN_FREE_GB:-4}
RUN_MIN_DISK_FREE_GB=${RUN_MIN_DISK_FREE_GB:-100}
INPUTS=${INPUTS:-}
OUTPUTS=${OUTPUTS:-}
DELETE_OUTPUTS=${DELETE_OUTPUTS:-0}
WORKDIR=${WORKDIR:-}
LEDGER=${LEDGER:-$HERE/ledger.tsv}
RESULTS=${RESULTS:-$HERE/runs}
NOTE=${NOTE:-}
HOST=$(hostname -s)
RUN=$RESULTS/${NAME}_$HOST
STDOUT=${STDOUT:-$RUN/stdout.txt}
(( $# > 0 )) || { echo "usage: NAME=.. MEM_LIMIT_GB=.. THREADS=.. run.bash <command> [args...]" >&2; exit 125; }
die() { echo "run.bash ABORT ($NAME): $*" >&2; exit 125; }

# Deletes $1 only if it is a file (not a folder: no recursive delete) strictly inside WORK_ROOT
# (resolved, no symlink games); anything else is refused and reported. The only way this script removes anything.
safe_rm() {
    local p root
    root=$(realpath -m "$WORK_ROOT")
    p=$(realpath -m "$1")
    if [[ $root == "${BENCH_DIR:-/tmp/colored-sbwt-bench}"/* && $p == "$root"/?* && $p != *"/../"* && -f $p ]]; then
        rm -- "$p"
    else
        echo "safe_rm: refusing to delete '$1' (resolved '$p'): not a file inside $root" >&2
        return 1
    fi
}

# MemAvailable plus the part of the ZFS ARC above c_min (see above); kernel_available_kb without the ARC
mem_available_kb() {
    awk '/^MemAvailable:/ { print $2 }' /proc/meminfo |
        awk -v arc=/proc/spl/kstat/zfs/arcstats '{ a = $1 }
            END { while ((getline l < arc) > 0) { split(l, f, " "); if (f[1] == "size") s = f[3]; if (f[1] == "c_min") m = f[3] }
                  print a + (s > m ? int((s - m) / 1024) : 0) }'
}
kernel_available_kb() { awk '/^MemAvailable:/ { print $2 }' /proc/meminfo; }
arc_kb() { awk '$1 == "size" { print int($3 / 1024) }' /proc/spl/kstat/zfs/arcstats 2> /dev/null; }
memory_pressure() { awk '/^full/ { split($2, f, "="); print int(f[2]) }' /proc/pressure/memory; }
# Busy and total jiffies of the whole server (busy = all but idle and iowait)
cpu_jiffies() { awk '/^cpu / { t = 0; for (i = 2; i <= 9; i++) t += $i; print t - $5 - $6, t }' /proc/stat; }
# Idle logical CPUs on average over $1 seconds
idle_cpus() {
    local b0 t0 b1 t1
    read -r b0 t0 < <(cpu_jiffies); sleep "$1"; read -r b1 t1 < <(cpu_jiffies)
    awk -v b=$(( b1 - b0 )) -v t=$(( t1 - t0 )) -v n="$(nproc)" 'BEGIN { printf "%.1f", n * (1 - b / t) }'
}
# Copies stdin, each line prefixed with the epoch in microseconds (as ts '%.s', absent here); its
# process is named stderr_ts so that the sampler can tell it from the command
stderr_ts() { perl -MTime::HiRes=time -ne 'BEGIN { $0 = "stderr_ts"; $| = 1 } printf "%.6f %s", time, $_'; }
others_top() { ps -eo user:32,pid,rss,pcpu,etime,args --sort=-rss --no-headers | awk -v me="$(id -un)" '$1 != me' | head -"$1"; }
disk_free_gb() { df --output=avail -B1G "$WORK_ROOT" | tail -1 | tr -d ' '; }
work_used_gb() { du -s -B1G "$WORK_ROOT" 2> /dev/null | cut -f1; }

# The files of INPUTS, one per line (folders: the files in them; @<file>: the paths in its first column)
input_files() {
    local i
    for i in $INPUTS; do
        if [[ $i == @* ]]; then cut -f1 "${i#@}"
        elif [[ -d $i ]]; then find "$i" -type f
        else echo "$i"; fi
    done
}

[[ -e $RUN ]] && die "$RUN exists already (delete it, or its ledger line, to run again)"
missing=$(input_files | while read -r f; do [[ -r $f ]] || echo "$f"; done | head -3)
[[ -z $missing ]] || die "input not readable: $missing"
mkdir -p "$WORK_ROOT" "$RESULTS"
[[ -n $WORKDIR ]] && { [[ $(realpath -m "$WORKDIR") == "$(realpath -m "$WORK_ROOT")"/?* ]] || die "WORKDIR $WORKDIR is not inside $WORK_ROOT"; }
for f in $OUTPUTS; do
    [[ $DELETE_OUTPUTS == 1 && $(realpath -m "$f") != "$(realpath -m "$WORK_ROOT")"/?* ]] && die "DELETE_OUTPUTS=1 but output $f is not inside $WORK_ROOT"
done

# ---- Wait for the lock, then for a quiet server
exec 9> "$LOCK"
flock -n 9 || { echo "run.bash ($NAME): waiting for the lock $LOCK since $(date -Is)" >&2; flock 9; }
(( $(disk_free_gb) >= MIN_DISK_FREE_GB )) || die "only $(disk_free_gb) GiB free on the disk of $WORK_ROOT, need MIN_DISK_FREE_GB=$MIN_DISK_FREE_GB"
(( $(work_used_gb) <= DISK_BUDGET_GB )) || die "$WORK_ROOT holds $(work_used_gb) GiB, over DISK_BUDGET_GB=$DISK_BUDGET_GB"
waited=0
while :; do
    avail_gb=$(( $(mem_available_kb) / 1048576 ))
    idle=$(idle_cpus 30)
    (( avail_gb >= MEM_LIMIT_GB + MIN_FREE_GB )) && awk -v i="$idle" -v m="$MIN_IDLE_CPUS" 'BEGIN { exit !(i >= m) }' && break
    (( WAIT_S > 0 && waited >= WAIT_S )) && die "server not quiet after ${waited} s: ${avail_gb} GiB available (need $(( MEM_LIMIT_GB + MIN_FREE_GB ))), $idle idle CPUs (need $MIN_IDLE_CPUS)"
    (( waited % 600 == 0 )) && echo "run.bash ($NAME): waiting for a quiet server: ${avail_gb} GiB available (need $(( MEM_LIMIT_GB + MIN_FREE_GB ))), $idle idle CPUs (need $MIN_IDLE_CPUS), $(date -Is)" >&2
    sleep 30; waited=$(( waited + 60 ))   # + the 30 s of idle_cpus
done

# ---- Memory cap, checked first on a probe and again inside the cgroup of the command
mem_max=$(( MEM_LIMIT_GB << 30 ))
cap=(systemd-run --user --scope --quiet -p MemoryMax="$mem_max" -p MemorySwapMax=0)
probe=$("${cap[@]}" sh -c 'cat "/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)/memory.max"' 2>&1)
[[ $probe == "$mem_max" ]] || die "no memory cap: systemd-run --user --scope gave '$probe', expected $mem_max"
in_cap='cg=$(cut -d: -f3 /proc/self/cgroup); max=$(cat "/sys/fs/cgroup$cg/memory.max")
echo "cgroup $cg memory.max=$max" >&2; [ "$max" = "$0" ] || exit 99; exec "$@"'

mkdir -p "$RUN"
{
    echo "name=$NAME host=$HOST start=$(date -Is) threads=$THREADS cap=${MEM_LIMIT_GB}G min_free=${MIN_FREE_GB}G min_idle_cpus=$MIN_IDLE_CPUS"
    echo "at start: ${avail_gb} GiB available, $idle idle CPUs, load $(cut -d' ' -f1-3 /proc/loadavg), disk free $(disk_free_gb) GiB, waited ${waited} s"
    echo "outputs=$OUTPUTS workdir=$WORKDIR delete_outputs=$DELETE_OUTPUTS note=$NOTE"
    printf 'command:'; printf ' %q' "$@"; echo
} > "$RUN/cmd.txt"

# ---- Inputs out of the page cache (cold start)
if [[ -n $INPUTS ]]; then
    input_files | tr '\n' '\0' | xargs -0 -r sync 9>&- || die "sync of the inputs failed"
    input_files | while read -r f; do dd if="$f" iflag=nocache count=0 status=none || exit 1; done 9>&- || die "dropping the inputs from the page cache failed"
    # (xargs may split the files over several calls: one line each, summed)
    read -r cached_in bytes_in < <(input_files | tr '\n' '\0' | xargs -0 -r python3 -I "$HERE/cached_bytes.py" | awk '{ c += $1; t += $2 } END { print c + 0, t + 0 }')
    printf 'files=%s bytes=%s still_cached=%s inputs=%s\n' "$(input_files | wc -l)" "$bytes_in" "$cached_in" "$INPUTS" > "$RUN/cache.txt"
    echo "cold start: $(cat "$RUN/cache.txt")" >> "$RUN/cmd.txt"
fi

# server.tsv, workdir.tsv, others.tsv and the memory watchdog, every 2 s while $1 runs (see above)
watchdog() {
    local pid=$1 cg reason n=0 n_low=0 b t ev pk
    printf 'epoch\tcgroup_kb\tanon_kb\tfile_kb\tcgroup_cpu_usec\tserver_busy_jiffies\tserver_total_jiffies\tavailable_kb\tpressure_full_avg10\tkernel_available_kb\tarc_kb\n' > "$RUN/server.tsv"
    printf 'epoch\tworkdir_kb\n' > "$RUN/workdir.tsv"
    printf 'epoch\tuser\tpid\trss_kb\tcpu_percent\telapsed\tcommand\n' > "$RUN/others.tsv"
    while [[ -d /proc/$pid ]]; do
        (( n % 15 == 0 )) && others_top 5 | awk -v t="$(date +%s)" -v OFS='\t' '{ cmd = $6; for (i = 7; i <= NF; i++) cmd = cmd " " $i; print t, $1, $2, $3, $4, $5, substr(cmd, 1, 200) }' >> "$RUN/others.tsv"
        [[ -n $WORKDIR ]] && printf '%s\t%s\n' "$(date +%s)" "$(du -s -k "$WORKDIR" 2> /dev/null | cut -f1)" >> "$RUN/workdir.tsv"
        n=$(( n + 1 ))
        cg=/sys/fs/cgroup$(cut -d: -f3 "/proc/$pid/cgroup" 2> /dev/null)
        if [[ $cg == */run-*.scope ]]; then
            # (only whole: at the end of the run the cgroup can disappear between the two reads)
            ev=$(cat "$cg/memory.events" 2> /dev/null) && pk=$(cat "$cg/memory.peak" 2> /dev/null) &&
                printf '%s\npeak %s\n' "$ev" "$pk" > "$RUN/cgroup.tmp" && mv "$RUN/cgroup.tmp" "$RUN/cgroup.txt"
            read -r b t < <(cpu_jiffies)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date +%s.%N)" "$(( $(cat "$cg/memory.current" 2> /dev/null || echo 0) / 1024 ))" \
                "$(awk '$1 == "anon" { print int($2 / 1024) }' "$cg/memory.stat" 2> /dev/null)" \
                "$(awk '$1 == "file" { print int($2 / 1024) }' "$cg/memory.stat" 2> /dev/null)" \
                "$(awk '$1 == "usage_usec" { print $2 }' "$cg/cpu.stat" 2> /dev/null)" "$b" "$t" \
                "$(mem_available_kb)" "$(awk '/^full/ { split($2, f, "="); print f[2] }' /proc/pressure/memory)" \
                "$(kernel_available_kb)" "$(arc_kb)" >> "$RUN/server.tsv"
        fi
        (( $(kernel_available_kb) < KERNEL_MIN_FREE_GB * 1048576 )) && n_low=$(( n_low + 1 )) || n_low=0
        reason=""
        if (( $(mem_available_kb) < MIN_FREE_GB * 1048576 )); then
            reason="server available memory under MIN_FREE_GB=$MIN_FREE_GB"
        elif (( n_low >= 5 )); then
            reason="kernel available memory under KERNEL_MIN_FREE_GB=$KERNEL_MIN_FREE_GB for 10 s (the ARC does not shrink fast enough)"
        elif (( $(memory_pressure) > MAX_PRESSURE && $(kernel_available_kb) < PRESSURE_FREE_GB * 1048576 )); then
            reason="memory pressure (full avg10) $(memory_pressure) % over MAX_PRESSURE=$MAX_PRESSURE with $(( $(kernel_available_kb) / 1048576 )) GiB available"
        elif (( $(disk_free_gb) < RUN_MIN_DISK_FREE_GB )); then
            reason="disk of $WORK_ROOT at $(disk_free_gb) GiB free, under RUN_MIN_DISK_FREE_GB=$RUN_MIN_DISK_FREE_GB"
        fi
        if [[ -n $reason ]]; then
            { echo "watchdog: $reason at $(date -Is)"
              others_top 3 | awk '{ cmd = $6; for (i = 7; i <= NF; i++) cmd = cmd " " $i; printf "  other: %s pid %s, %.1f GiB, %s%% CPU: %s\n", $1, $2, $3 / 1048576, $4, substr(cmd, 1, 120) }'
            } | tee "$RUN/killed.txt" >&2
            kill -TERM "$pid" 2> /dev/null; sleep 2; kill -KILL "$pid" 2> /dev/null
            break
        fi
        sleep 2
    done
}

pid=""
trap '[[ -n $pid ]] && kill -TERM "$pid" 2> /dev/null; echo "interrupted at $(date -Is)" >> "$RUN/killed.txt"; exit 130' INT TERM HUP

# Every command up to /usr/bin/time execs the next one, so $! is the pid of time; the command is its child
"${cap[@]}" sh -c "$in_cap" "$mem_max" choom -n 1000 -- nice -n 19 ionice -c 2 -n 7 \
    /usr/bin/time -v -o "$RUN/time.txt" "$@" > "$STDOUT" 2> >(stderr_ts > "$RUN/stderr.log") 9>&- &
pid=$!
cmd_pid=""
# (time has two children: the command, and stderr_ts of the process substitution of stderr)
for _ in $(seq 100); do
    cmd_pid=$(ps -o pid=,comm= --ppid "$pid" | awk '$2 != "stderr_ts" { print $1; exit }')
    [[ -n $cmd_pid ]] && break; [[ -d /proc/$pid ]] || break; sleep 0.05
done
[[ -n $cmd_pid ]] && "$HERE/sample_proc.bash" "$cmd_pid" "$RUN/samples.tsv" "$SAMPLE_S" 9>&- &
sampler=$!
watchdog "$pid" 9>&- &
dog=$!
wait "$pid"; status=$?
pid=""
wait "$sampler" 2> /dev/null
kill "$dog" 2> /dev/null; wait "$dog" 2> /dev/null
sleep 0.5   # let stderr_ts flush stderr.log
grep -m1 '^[0-9.]* cgroup ' "$RUN/stderr.log" >> "$RUN/cmd.txt"
(( status == 99 )) && die "command not started: $(head -1 "$RUN/stderr.log"), expected memory.max=$mem_max"
[[ -f $RUN/killed.txt ]] && status=137

# ---- Results: /usr/bin/time, cgroup peak, CPU of the other users, outputs
tv() { awk -F': ' -v k="$1" 'index($0, k) { print $2; exit }' "$RUN/time.txt" 2> /dev/null; }
wall=$(tv 'Elapsed (wall clock) time' | awk -F: '{ s = 0; for (i = 1; i <= NF; i++) s = s * 60 + $i; printf "%.2f", s }')
user=$(tv 'User time (seconds)'); sys=$(tv 'System time (seconds)'); maxrss=$(tv 'Maximum resident set size (kbytes)')
cg_peak=$(awk '$1 == "peak" { print $2 }' "$RUN/cgroup.txt" 2> /dev/null)
oom=$(awk '$1 == "oom_kill" { print $2 }' "$RUN/cgroup.txt" 2> /dev/null)
# Mean logical CPUs busy with other processes during the run: server busy minus this cgroup
others_cpu=$(awk -F'\t' -v hz="$(getconf CLK_TCK)" 'NR == 2 { b0 = $6; c0 = $5; t0 = $1 } NR > 2 { b = $6; c = $5; t = $1 }
    END { if (t > t0) printf "%.1f", ((b - b0) / hz - (c - c0) / 1e6) / (t - t0); else print "NA" }' "$RUN/server.tsv" 2> /dev/null)
work_peak_kb=$(awk -F'\t' 'NR > 1 && $2 > m { m = $2 } END { print m + 0 }' "$RUN/workdir.tsv" 2> /dev/null)
sizes="" md5s=""
for f in $OUTPUTS; do
    if [[ -f $f ]]; then
        sizes+="$(stat -c %s "$f"),"; md5s+="$(md5sum < "$f" | cut -d' ' -f1),"
    else
        sizes+="NA,"; md5s+="NA,"
    fi
done
killed=$(head -1 "$RUN/killed.txt" 2> /dev/null)
header='name\thost\tstart\texit\twall_s\tuser_s\tsys_s\tmax_rss_kb\tcgroup_peak_bytes\toom_kill\tthreads\tcap_gb\tothers_cpu_mean\tworkdir_peak_kb\toutput_bytes\toutput_md5\tkilled\tnote'
line=$(printf '%s\t' "$NAME" "$HOST" "$(awk '{ sub(/^start=/, "", $3); print $3 }' "$RUN/cmd.txt" | head -1)" "$status" "$wall" "$user" "$sys" "$maxrss" \
    "${cg_peak:-NA}" "${oom:-NA}" "$THREADS" "$MEM_LIMIT_GB" "$others_cpu" "$work_peak_kb" "${sizes%,}" "${md5s%,}" "$killed" | sed 's/\t$//')
line+=$'\t'"$NOTE"
{ printf "$header\n"; echo "$line"; } > "$RUN/result.txt"
(
    flock 8
    [[ -s $LEDGER ]] || printf "$header\n" > "$LEDGER"
    echo "$line" >> "$LEDGER"
) 8>> "$LEDGER.lock"
printf "run.bash ($NAME): exit %s, wall %s s, CPU %s+%s s, peak RSS %.2f GiB, others %s CPUs\n" "$status" "$wall" "$user" "$sys" \
    "$(awk -v k="${maxrss:-0}" 'BEGIN { print k / 1048576 }')" "$others_cpu" >&2

if [[ $DELETE_OUTPUTS == 1 ]]; then for f in $OUTPUTS; do [[ -e $f ]] && safe_rm "$f"; done; fi
if [[ -n $WORKDIR && -d $WORKDIR ]]; then
    for f in "$WORKDIR"/* "$WORKDIR"/.[!.]*; do [[ -e $f ]] && safe_rm "$f"; done
fi
exit "$status"
