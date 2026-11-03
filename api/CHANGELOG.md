# Changelog

## Unreleased

### New features

* Added `merge_without_cleanup`, which merges two SBWTs like `merge`, but keeps the dummy nodes that become redundant in the merge. The result has the same k-mers but possibly more dummy nodes, and its colex positions correspond exactly to the positions of the `MergeInterleaving`. `merge` is unchanged.
* Added `merge_with_removed_positions`, which merges like `merge`, but also returns which positions of the `MergeInterleaving` were removed.
* Added the flag `--keep-redundant-dummies` to the CLI command `merge`, which uses `merge_without_cleanup`.
* Added `intersect_with_removed_positions`, which intersects like `intersect`, but also returns the interleaving the result consists of (a `ThreeWayInterleaving`, whose third set of nodes is the dummy nodes of the repaired dummy chains) and which of its positions are not in the result. `ThreeWayInterleaving` is now public.
* Added `difference_with_removed_positions`, which does for `difference` what `intersect_with_removed_positions` does for `intersect`.
* Added `filter` and `filter_with_removed_positions`, which keep the k-mers of an index that are marked in a bit vector. The result is the same as an index built from scratch from those k-mers. No LCS array is needed.
* Added `result_k_minus_1_marks`, which gives the (k-1)-suffix groups of the result of a set operation from the merge plan that the `_with_removed_positions` variants return, and `Dbg::new_with_k_minus_1_marks`, which builds a DBG with given groups. Together they build the DBG of such a result without an LCS array or a column inversion.

### Changes

* `intersect`, `difference` and `filter` repair the dummy chains of their result in the given `MergeInterleaving`, instead of building an auxiliary index of the affected k-mers and a three-way interleaving of the inputs and that index from scratch. This is faster and needs less memory, with the same result. Their `optimize_peak_ram` parameter no longer has an effect.
* `intersect`, `difference` and `filter` remove the dummy chains that lead to no k-mer of the result in parallel, with one bit per node instead of a list of the removed nodes. This is much faster and needs much less memory on inputs with many short sequences (e.g. short reads).
* `merge` now also removes redundant dummy nodes that were already in its inputs, e.g. in the output of `merge_without_cleanup`, so the result is always the same as an index built from scratch. Finding the redundant dummy nodes is also faster.
* The removal of redundant dummy nodes in `merge` uses less memory and is faster, especially on inputs with many dummy nodes (e.g. short reads). 
* `Dbg::new` without an LCS array is faster: the column inversion that marks the (k-1)-suffix groups now compares the columns in parallel. `filter` without the groups gets the same speedup.
* `Dbg::walk_unitig_from`, and with it `Dbg::iter_unitigs_with_callback`, is faster: one rank per node instead of a rank and a select. Its `out_labels_buf` parameter is no longer used.

### Bugfixes

* `MergeInterleaving::new`, `merge`, `intersect` and `difference` no longer panic when there are more threads than positions to split between them, which happened on small inputs.

## v0.6.3

### New features

* Added a new subset rank structure: SubsetCorrectionSets. It uses roughly half the space compared to SubsetMatrix, but is 4-5 times slower.
* Added an option to stream the suffix array from disk (instead of keeping it in memory) to `BuildByBoundedSuffixSort` and the libsais-based builder, reducing peak memory usage during construction.

### Bugfixes

* `SbwtIndex::from_subset_seq` now builds rank support automatically if it is missing, instead of assuming it is already present.
* Fixed a buffer-length bug in `par_build_lengths`.
* Fixed `get_C_array` to also work on subset sequences without rank support.
* Various edge case fixes in the rank/select data structures (`pred8_vs1`, `quad_rank`, `pino`, `subsetseq`), found via added edge-case tests.

### Performance 

* Improved the peak space usage during in-memory construction.

## v0.6.2

Removed some unnecessary implementation details from the docs.

## v0.6.1

The BuildByLibsais construction algorithm now runs fully in linear time at the cost of using O(n log n + n log k) extra bits. BuildByBoundedSuffixSort uses the same functions, so it also pays the same memory cost.

## v0.6.0

### New features

* A variable-order DBG interface, with optional k-mer counts.
* Two new construction algorithms: one based on bounded-context suffix sorting, another based on suffix sorting via libsais (behind the feature gate "libsais"). The builder API has been reworked to accommodate these, and it is **not backwards compatible**.

### Bugfixes

* Fixed a bug in the set difference operation which lead to a corrupt index.
* Fixed a bug where redundant dummy nodes were sometimes left behind by set operations.
* Fixed a crash when running set operations on very small inputs.

### Performance

* Disk-based construction: better batching, and introduced parallelism for later stages of the construction pipeline.

## v0.5.0

- Added support for set operations (union, intersection, difference).
  - API: new `set_operations` module (`src/set_operations`) exposing the
    core functions `merge`, `intersect` and `difference`, together with
    helper types such as `MergeInterleaving`.
  - Union is based on the existing merge code. It now also cleans up redundant dummy nodes.
- Fixed bug in `push_labels_forward_compact` which occurred in some cases when some nucleotide does not occur in the SBWT at all.

## v0.4.3

* Update to jseqio v0.1.7. This required adding the Sync trait bound to `run_from_fasta` and `run_from_fastq` in the builder.

## v0.4.2

### Performance

* Significantly speed up `push_labels_forward`, speeding up e.g. the SBWT merge algorithm.
* Better parallelism in unitig iteration and exporting.

### API features

* Expose a function for a single matching statistics update step.
* Add a matching statistics iterator with a length bound below k.
* Make functions related to unitigs public.
* Make the dummy marking algorithm public and move it to SbwtIndex.
* Provide access to DBG dummy marks.
* Eliminate most dependencies on `simple_sds_sbwt` from the public API. Some dependencies remain at `streaming_index.rs` and `sdsl_compatibility.rs`.

### Bugfixes

* Make the DBG implementation work even if there are redundant dummies.
* Fix empty ranges in parallel segmentation.
* Fix past-the-end bits in bitvec.

## v0.4.1

New algorithms and performance improvements. This is a major update, but still backward compatible.

* Add an algorithm to merge two SBWTs.
* Improvements in the parallelism and performance of the in-memory construction algorithm.
* In-memory construction is no longer behind a feature gate.
* Rename BitPackedKmerSorting to BitPackedKmerSortingDisk. The old type name works as an alias but is marked as deprecated.
* Add a streaming matching statistics iterator.

## v0.4.0

- Added `iter_unitigs_with_callback`. This can break existing code because now `parallel_export_unitigs` also takes the
  number of threads as an input parameter, whereas previously it used all available threads.

## v0.3.11

- Fix SbwtIndex documentation not showing up.

## v0.3.10

- Fix `push_labels_forward` behavior at root: now fills in dollars, as documentation says.

## v0.3.9

- Make `build_last_column` public.

## v0.3.8

- Add disk usage reporting during construction.
- Make `push_labels_forward` public.

## v0.3.7

- Parallelise the subroutine `push_labels_forward`. This speeds up the k-mer dump, LCS array construction and DBG construction.
- Add a CLI command `build-lcs` to build the LCS array from the SBWT on the command line.
- Make CLI commands `lookup`, `matching-statistics`, `dump-kmers` and `dump-unitigs` all print to stdout if an output file is not given.

## v0.3.6

- Fixed a corner case when importing an index with an empty precalc table from the C++ format. This caused a crash at query time, which is fixed in this patch.


## v0.3.5

- Add `load_from_cpp_plain_matrix_format` to load the index from the C++ API
  or CLI format: https://github.com/algbio/SBWT.

## v0.3.4

- Add an optional in-memory construction implementation by Tommi Mäklin: https://github.com/jnalanko/sbwt-rs-cli/pull/4. This can be built with by enabling the feature `bpks-mem`.
- Fixes to make the crate work with WebAssembly.
- Parallel C-array construction.

## v0.3.3

- Documentation: Add an example using Needletail for FASTX parsing.
- Change usize to u32 in `IS_DNA` for platform-independence.

## v0.3.2

- Fixes to documentation.

## v0.3.1

- Added function `from_subset_seq` to SbwtIndex.

## v0.3.0

- Added serialization magic string to the start of the serialization format. This breaks compatibility with files generated using the previous versions.
- Made some functions generic over the particular implementation of SubsetSeq in preparation for the future when we will have multiple implementations.
- Added to the documentation a link to the Github repository of the CLI.

## v0.2.1

- Added CHANGELOG.md
- Updated jseqio to version 0.1.3. This makes it so that lower-case characters are made upper-case during parsing. Previously k-mers containing lower-case characters were ignored.
