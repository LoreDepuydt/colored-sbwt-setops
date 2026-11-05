# colored-sbwt-setops

Edit a colored k-mer index without rebuilding it and without the original genomes: remove or extract
genomes or species, or compute the union, intersection or difference of two indexes. The result is a
complete colored index, ready to query, with the same k-mers and color sets as the index built from
scratch from the genomes it contains.

The indexes are those of [Themisto 2](https://github.com/jnalanko/themisto2),
and are based on the spectral Burrows-Wheeler transform (SBWT). This repository is Themisto 2 and the
[sbwt](https://github.com/jnalanko/sbwt-rs-cli) library with these operations added. It has all commands
of Themisto 2, so you can build, edit and query an index with one tool.

## When to use it

- Genomes turn out to be contaminated or are retracted: remove them from the index.
- You need an index of one species or clade: extract it from a large index.
- New genomes arrive: build an index of them and merge it into the existing one.
- You want to compare two collections: find the k-mers that they share, or that only one of them has.

Editing saves the most when the result is a large part of the index. Removing one genome from an index
of 3,682 *E. coli* genomes took 54 seconds instead of 13 minutes for a rebuild, and adding one took
about a minute. Extracting *M. tuberculosis* (180,094 genomes) from the AllTheBacteria index of 12,733
species and 192 billion k-mers took less than four hours.

## Installation

Requires Rust 1.89 or later, which you can install with [rustup](https://rust-lang.org/tools/install/)
without root privileges.

```bash
git clone https://github.com/LoreDepuydt/colored-sbwt-setops
cd colored-sbwt-setops/colors
cargo build --release
export PATH="$PWD/target/release:$PATH"
```

The binary is `colors/target/release/themisto2`. It has this name because this is a development version
of Themisto 2, into which these operations will be merged. It is not an official release.

## Quick start

Run this in `colors/`, on the small genomes in `example/`. Their sequences have 25 to 55 bases, hence the
small k: a sequence shorter than k has no k-mers, so it adds nothing to the index. For real genomes use
k = 31.

Build an index from a file with one line per color, `<FASTA file> TAB <color name>`, and remove a color:

```bash
printf 'example/C1.fna\tC1\nexample/C2.fna\tC2\nexample/C3.fna\tC3\n' > colors.tsv
themisto2 build -k 5 --file-colors colors.tsv -o all.thm2 --temp-dir tmp -t 4
themisto2 subset --index all.thm2 --colors C1 --complement -o no_c1.thm2 -t 4
```

Query the result as any index, here for the colors that have at least 70% of the k-mers of each query:

```bash
themisto2 threshold-pseudoalign --threshold 0.7 --denominator all --min-hits 1 -i no_c1.thm2 -q example/C1.fna -t 4
```

```
{"name": "C1.1", "colors": [0]}
{"name": "C1.2", "colors": [0,1]}
```

Colors are numbered per index: color 0 of `no_c1.thm2` is C2. To see the names, run `themisto2
dump-color-names -i no_c1.thm2`, or give the output of the query to `themisto2 report`.

## Commands

Every command prints its options with `--help`.

```bash
# Extract genomes: the k-mers of C1 and C2, with only these colors
themisto2 subset --index all.thm2 --colors C1,C2 -o c12.thm2 -t 4
# ... with the names in a file, one per line
printf 'C1\nC2\n' > keep.txt
themisto2 subset --index all.thm2 --colors-file keep.txt -o c12.thm2 -t 4
# ... or the k-mers of C3, each with all colors of the index that have it
themisto2 subset --index all.thm2 --colors C3 --keep-all-colors -o c3_all.thm2 -t 4

# Union: merge any number of indexes, listed in a file with one index per line
themisto2 subset --index all.thm2 --colors C3 -o c3.thm2 -t 4
printf 'c12.thm2\nc3.thm2\n' > indexes.txt
themisto2 merge --index-file-list indexes.txt --temp-dir tmp -o merged.thm2 -t 4

# Intersection: the k-mers in both indexes, by default with their colors of either index
themisto2 intersect --index1 all.thm2 --index2 c3.thm2 --merge-shared-colors -o inter.thm2 -t 4
# ... or with only the colors that are in both
themisto2 intersect --index1 all.thm2 --index2 c3.thm2 --colors intersect -o inter_colors.thm2 -t 4

# Difference: the k-mers of the first index that are not in the second
themisto2 difference --index1 all.thm2 --index2 c3.thm2 -o diff.thm2 -t 4

# Inspect an index: statistics, color names, or the unitigs (GFA) with their color sets
themisto2 stats -i diff.thm2
themisto2 dump-color-names -i diff.thm2
themisto2 export -i diff.thm2 -o diff --gfa -t 4
```

Good to know:

- **Color names.** `subset` selects colors by name. If the list given to `build` has no names, the name
  of a color is its path as written in the list, for example `example/C1.fna`.
- **Colors with the same name.** `merge` and `intersect` keep the colors of the two inputs apart, so a
  name can occur twice in the result. To make them one color, pass `--merge-shared-colors`, which is
  what you want when the inputs come from the same collection. Duplicate names can not be merged later.
- **The same k.** All inputs of an operation must have the same k.
- **Output files.** Use a new name for `-o`: `build`, `merge`, `intersect` and `difference` clear the
  output when they start, so an output that is also an input is lost.
- **Memory and disk.** The set operations and queries load the whole index into memory.
  `--low-ram-mode` uses slower algorithms that need less memory; the result is the same. `export` writes
  text, which can be far larger than the index.
- **Export and import.** `import` numbers the colors: the color names are not kept.
- **The LCS array.** Queries use it to move from one k-mer to the next in about one step per base. An
  index without one (`--no-lcs`) can still be queried, but each k-mer is searched from scratch: about
  1.5 to 1.75 times slower in a small test, and probably more on large indexes.

## Query modes

For each query sequence, both modes write one JSON line with the name of the sequence and its colors.

- `intersection-pseudoalign` reports the colors of *every* k-mer of the query that is in the index.
- `threshold-pseudoalign` reports the colors of at least a fraction of the k-mers of the query, given by
  `--threshold`. This is more robust to sequencing errors, small variants and incomplete references, but
  less sensitive to small strain-level variation.

## Large indexes

For more than a few gigabytes of repetitive genomes, deduplicate the k-mers first. The final index is the
same, and the construction can be orders of magnitude faster. This needs
[GGCAT](https://github.com/algbio/ggcat) (BCALM2 and Cuttlefish 3 also work) and the `sbwt` tool of this
repository: run `cargo build --release` in `sbwt-rs-cli/`.

```bash
# Unitigs of all genomes together, without colors. GGCAT needs a list with only the file names.
cut -f1 colors.tsv > genomes.txt
ggcat build -l genomes.txt -o unitigs.fna -s 1 -k 31 -t temp -j 32 -m 64 -p
# The SBWT and LCS array of the unitigs: unitigs.sbwt and unitigs.lcs
sbwt build -i unitigs.fna -o unitigs -k 31 -r -l -t 32 -m 64 -v --temp-dir temp
# The colors, from the original genomes
themisto2 build --file-colors colors.tsv -s unitigs.sbwt -l unitigs.lcs -o index.thm2 -k 31 -t 32
```

These commands use 32 threads and 64 GiB of memory (`-m`). In `ggcat`, `-s 1` is very important: without
it, GGCAT discards the k-mers that occur once. The temporary folder of `sbwt` should be on a fast disk with
up to 8 bytes per distinct k-mer. Keep the LCS array (`-l` in `sbwt build` and `themisto2 build`), to make
queries faster.

The SBWT must be built with `-r`, from exactly the genomes of the color list: the k-mers of a genome that
is only in the SBWT silently get no colors.

## Repository layout

| Folder | What |
|---|---|
| `colors/` | the colored index and the operations on it, the command-line tool `themisto2` (from [themisto2](https://github.com/jnalanko/themisto2)) |
| `sbwt-rs-cli/` | the SBWT library and its command-line tool, with the operations on SBWTs (from [sbwt-rs-cli](https://github.com/jnalanko/sbwt-rs-cli)) |
| `benchmarks/` | the scripts of the experiments of the paper |

Both tools will be merged into their upstream repositories in time.

## Citation

If you use this software, please cite:

> Lore Depuydt, Jarno N. Alanko, Camille Marchet, and Simon J. Puglisi. *Extracting and Updating Compact
> Colored k-mer Indexes.* 2026.

A preprint will follow.

The construction algorithm of Themisto 2 is described in:

> Jarno N. Alanko and Simon J. Puglisi. [*Construction of Distinct k-mer Color Sets via Set
> Fingerprinting.*](https://doi.org/10.4230/LIPIcs.WABI.2026.6) WABI 2026, LIPIcs 390, 6:1–6:17.

## Reproducing the paper

The scripts, and how to run them, are in [benchmarks/README.md](benchmarks/README.md). The code of the
paper is pinned by two tags:

| Tag | Used for |
|---|---|
| `recomb27` | all results of the paper |
| `recomb27-threeway` | the comparison with the three-way merge of our earlier work, which this tag puts back in place of insertion |

## License

MIT, see [LICENSE](LICENSE). `sbwt-rs-cli/` keeps its own MIT license, in
[sbwt-rs-cli/LICENSE.txt](sbwt-rs-cli/LICENSE.txt).
