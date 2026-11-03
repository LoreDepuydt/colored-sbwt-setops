# colored-sbwt-setops

Set operations on colored k-mer indexes based on the spectral Burrows-Wheeler transform (SBWT):
extracting genomes or species from an index, and the union, intersection and difference of two
indexes. Each operation works on the indexes alone, without the genomes, and gives a complete colored
index, the same as the index built from scratch from the genomes it contains.

This is the implementation of the paper *Extracting and Updating Compact Colored k-mer Indexes*.

## Contents

| Folder | What |
|---|---|
| `sbwt-rs-cli/` | the SBWT library and its command-line tool, with the operations on SBWTs (from [sbwt-rs-cli](https://github.com/jnalanko/sbwt-rs-cli)) |
| `colors/` | the colored index and the operations on it, the command-line tool `themisto2` (from [themisto2](https://github.com/jnalanko/themisto2)) |
| `benchmarks/` | the scripts of the experiments of the paper, see [benchmarks/README.md](benchmarks/README.md) |

Both will be merged into their upstream repositories.

## Build

Requires Rust 1.77 or later.

```bash
cd colors
cargo build --release
```

The binary is `colors/target/release/themisto2`.

## Usage

Every command prints its options with `--help`. An example on the three small genomes in
`colors/example/` (run in `colors/`, with `themisto2` on the `PATH`):

```bash
# Build an index, with one color per genome: a file of <FASTA file> TAB <color name>
printf 'example/C1.fna\tC1\nexample/C2.fna\tC2\nexample/C3.fna\tC3\n' > colors.tsv
themisto2 build -k 31 --file-colors colors.tsv -o all.thm2 --temp-dir tmp -t 4

# Extract genomes: the k-mers of C1 and C2, with only these colors, as the index built from them
themisto2 subset --index all.thm2 --colors C1,C2 -o c12.thm2 -t 4
# ... or of all colors but C1 (--complement), or with the names in a file (--colors-file)
themisto2 subset --index all.thm2 --colors C1 --complement -o not_c1.thm2 -t 4
# ... or the k-mers of C3, each with all colors of the index that have it
themisto2 subset --index all.thm2 --colors C3 --keep-all-colors -o c3_all.thm2 -t 4

# Union: merge indexes (a file with one index per line)
themisto2 subset --index all.thm2 --colors C3 -o c3.thm2 -t 4
printf 'c12.thm2\nc3.thm2\n' > indexes.txt
themisto2 merge --index-file-list indexes.txt --temp-dir tmp -o merged.thm2 -t 4

# Intersection and difference of two indexes
themisto2 intersect --index1 all.thm2 --index2 c3.thm2 -o inter.thm2 -t 4
themisto2 intersect --index1 all.thm2 --index2 c3.thm2 --colors intersect -o inter_colors.thm2 -t 4
themisto2 difference --index1 all.thm2 --index2 c3.thm2 -o diff.thm2 -t 4

# Inspect a result: statistics, color names, or the unitigs (GFA) with their color sets
themisto2 stats -i diff.thm2
themisto2 dump-color-names -i diff.thm2
themisto2 export -i diff.thm2 -o diff --gfa -t 4
```

Options that the operations share:

- `--low-ram-mode`: slower algorithms that need less memory; the result is the same.
- `--merge-shared-colors` (merge, and intersect with `--colors union`): colors with the same name in
  the inputs are one color.
- `--no-lcs`: for benchmarking only, as in the paper; the result cannot be queried.

The results can be queried as any index of `themisto2` (`threshold-pseudoalign`,
`intersection-pseudoalign`).

## Versions of the paper

| Tag | Used for |
|---|---|
| `recomb27` | all results of the paper |
| `recomb27-threeway` | the comparison with the three-way merge of our earlier work, which this tag puts back in place of insertion |

## License

MIT, see [LICENSE](LICENSE). `sbwt-rs-cli/` keeps its own MIT license, in
[sbwt-rs-cli/LICENSE.txt](sbwt-rs-cli/LICENSE.txt).
