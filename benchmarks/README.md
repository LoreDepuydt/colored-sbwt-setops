# Benchmarks

The scripts of the experiments of the paper *Extracting and Updating Compact Colored k-mer Indexes*.

| Folder | Paper |
|---|---|
| `experiments/` | the benchmarks on *E. coli* and *S. pneumoniae*: Tables `exp-ecoli` and `exp-filter`, and the supplementary tables `supp-full`, `supp-steps` and `supp-pairs` |
| `atb/` | the case study on AllTheBacteria: its section, and the supplementary tables `supp-atb` and `atb` |

All experiments ran on a shared server with two Xeon Gold 6430 processors (64 cores, without
hyper-threading), 503 GiB of RAM and a local NVMe SSD, with 32 threads. Timings on other hardware will
differ; the results of the operations (the indexes) do not.

The scripts assume such a server: `run.bash` needs 200 GiB of free disk and waits for the memory of the
run plus 64 GiB to be available, and `download_pneumo.bash` needs 500 GiB of free disk on top of the download.

## Requirements

- Linux with cgroups v2 and a systemd user instance: every run is capped with
  `systemd-run --user --scope -p MemoryMax=...`.
- bash, Python 3 (standard library only), GNU coreutils, `/usr/bin/time`, `curl`, `md5sum`.
- [GGCAT](https://github.com/algbio/ggcat) 2.2.0, for the rebuilds.
- The binaries of this repository, in `experiments/bin/` (the setting `BIN` points elsewhere):

```bash
cd colors && cargo build --release --features extra-tools --bin themisto2 --bin compare_color_dumps && cd ..
cd sbwt-rs-cli && cargo build --release && cd ..
mkdir -p benchmarks/experiments/bin
cp colors/target/release/{themisto2,compare_color_dumps} sbwt-rs-cli/target/release/sbwt benchmarks/experiments/bin/
ln -s "$(command -v ggcat)" benchmarks/experiments/bin/ggcat
```

For the comparison with the three-way merge of our earlier work (`queue_oldfilter.bash`), also build
tag `recomb27-threeway` and copy its `themisto2` to `benchmarks/experiments/bin/oldfilter/themisto2`.

## Data

Everything large goes to one folder on a local disk, `BENCH_DIR` (default `/tmp/colored-sbwt-bench`):

| Path under `BENCH_DIR` | What | How |
|---|---|---|
| `data/ecoli/genomes/*.fna` | the 3,682 *E. coli* assemblies | from [Zenodo 10.5281/zenodo.6577997](https://doi.org/10.5281/zenodo.6577997) |
| `data/pneumo/genomes/*.fa.gz` | the *S. pneumoniae* assemblies of AllTheBacteria | `bash download_pneumo.bash <set>`, from the public bucket of AllTheBacteria, with md5 checks |
| `indexes/ecoli/scale_3682.thm2` | the index of all *E. coli* genomes | built by `queue_subset.bash` if missing |
| `indexes/pneumo/s_pneumoniae.thm2` | the index of the 121,034 *S. pneumoniae* genomes (43 GB) | see below |
| `work/` | outputs and temporary files of the runs | written by the scripts; `run.bash` deletes files only here |

The genome sets are in `experiments/sets/`: lists of color names, all cut from one random order of the
genomes (seed 1), made by `make_sets_ecoli.py` and `make_sets_pneumo.py`. `sets/pneumo/color_names.tsv`
holds the colors of `s_pneumoniae.thm2`, in order; `sets/pneumo/atb.tsv` the AllTheBacteria accession,
size and md5 of each genome.

To build `s_pneumoniae.thm2` from the assemblies, with the pipeline of the rebuilds (k = 31):

```bash
cd benchmarks/experiments
export BENCH_DIR=/tmp/colored-sbwt-bench
bash download_pneumo.bash all                                       # ~80 GB
G=$BENCH_DIR/data/pneumo/genomes T=$BENCH_DIR/work/pneumo_index && mkdir -p $T/tmp $BENCH_DIR/indexes/pneumo
awk -F'\t' -v g=$G '{ a = $2; sub(/.*\//, "", a); sub(/\.fa$/, "", a); print g "/" a ".fa.gz\t" $2 }' sets/pneumo/color_names.tsv > $T/fof.tsv
cut -f1 $T/fof.tsv > $T/genomes.txt
bin/ggcat build -l $T/genomes.txt -s 1 -k 31 -m 8 -p -j 32 -t $T/tmp -o $T/unitigs.fa
bin/sbwt build -i $T/unitigs.fa -o $T/pneumo -k 31 -t 32 -m 8 -r -v --temp-dir $T/tmp
bin/themisto2 build --file-colors $T/fof.tsv -s $T/pneumo.sbwt -o $BENCH_DIR/indexes/pneumo/s_pneumoniae.thm2 -k 31 -d 30 -t 32 --temp-dir $T/tmp --no-lcs
```

It has the colors of the index of the paper in the same order, and the same statistics (`themisto2
stats`: k-mers, colors, distinct color sets, unitigs); its extraction to the 1,000 genomes of
`sets/pneumo/keep_1000.txt` is the same as that of the index of the paper (`queue_rebuild.bash` checks
both).

## Running

Each `queue_*.bash` runs one group of experiments. Every run goes through `run.bash`, which waits for a
lock (one run at a time) and for a quiet server, drops the inputs from the page cache (every run starts
cold), caps the memory, samples the process, and writes one line to `ledger.tsv` (wall time, peak memory,
md5 of the outputs) and its files to `runs/<run>_<host>/`. A run that is in the ledger already is skipped,
so a queue can be stopped and started again. The header of each script lists its runs and settings.

```bash
cd benchmarks/experiments
export BENCH_DIR=/tmp/colored-sbwt-bench
bash download_pneumo.bash keep_12103      # the genomes of the S. pneumoniae rebuilds
bash pilot_ecoli.bash                     # optional: one run of each step, with checks
bash queue_subset.bash                    # extraction by colors vs rebuild; thread scaling, --low-ram-mode
bash queue_merge.bash                     # merge vs rebuild; builds the inputs of queue_setops.bash
bash queue_setops.bash                    # intersection and difference
bash queue_rebuild.bash                   # the full S. pneumoniae rebuilds
bash queue_parts.bash                     # the build time of the parts that are merged
bash queue_oldfilter.bash                 # the three-way merge of our earlier work against insertion
```

| Paper | Runs | Summary |
|---|---|---|
| Table `exp-ecoli`, `supp-full` | `queue_subset.bash`, `queue_merge.bash`, `queue_rebuild.bash`, `queue_parts.bash` | medians of `wall_s` and `max_rss_kb` in `ledger.tsv`; `python3 parts_summary.py` for the parts |
| Table `supp-steps` | `queue_subset.bash` | `python3 subset_log_stats.py` (time per step, from the logs) |
| Tables `exp-filter`, `supp-pairs` | `queue_oldfilter.bash`, with `queue_subset.bash` and `queue_setops.bash` | `python3 oldfilter_stats.py` |

The checks that the results equal the index built from scratch are part of the queues (`checks.tsv`);
`merge_check.bash`, `stream_check.bash` and `pneumo_check.bash` are the full checks on *S. pneumoniae*.
All runs use `--no-lcs`, as in the paper: the operations then find the (k-1)-suffix groups from the
SBWTs, and queries on the results are slower, as they search each k-mer from scratch.

## Case study (`atb/`)

`atb/mtb_case.bash` runs every step of the case study (its header lists them) through
`atb/capped_run.bash`, and `atb/m_direct.bash` the extraction of the *M. tuberculosis* color alone.
`atb/mtb2025/` builds the index of the *M. tuberculosis* genomes of releases 2024-08 and 2025-05 from
their assemblies (`make_lists.py`, then `download.bash`, then `MODE=all bash mtb2025.bash`).
`atb/gfastat/` is the tool for the unitig statistics (length, GC) of a GFA export
(`cargo build --release`). The inputs are settings, listed in the header of each script:

- the species index of AllTheBacteria, release 2024-08 (192 billion k-mers, 194 GB);
- the AllTheBacteria metadata (`ATB_META`): species calls, sylph and the file lists of releases 2024-08
  and 2025-05, and the GTDB r214 taxonomy;
- an index of the reference genomes H37Rv (NC_000962) and *M. bovis* AF2122/97 (NC_002945), one color,
  and one of the human genome CHM13 v2.0, one color: `themisto2 build -k 31 --file-colors <fof>`.
