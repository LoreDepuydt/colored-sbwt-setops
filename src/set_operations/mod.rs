//! Set operations on colored indexes: merging (union), intersecting and subtracting two indexes,
//! and taking the subset of an index for some of its colors.

pub(crate) mod colors;
pub(crate) mod key_kmers;
pub(crate) mod result_positions;
pub(crate) mod color_sets;
pub(crate) mod pipeline;
pub(crate) mod merge;
pub(crate) mod intersect;
pub(crate) mod difference;
pub(crate) mod subset;
pub(crate) mod dbg_cache;
#[cfg(test)]
mod test_utils;

use std::sync::atomic::{AtomicBool, Ordering};

use sbwt::{dbg::Dbg, LcsArray, SubsetMatrix};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

/// Whether the set operations use LCS arrays. Without (`--no-lcs`, for benchmarking only), their
/// inputs are loaded without their LCS arrays, every DBG computes its (k-1)-suffix groups from its
/// SBWT instead, and the result gets an empty LCS array, so that it cannot be queried.
static USE_LCS: AtomicBool = AtomicBool::new(true);

pub(crate) fn set_use_lcs(use_lcs: bool) {
    USE_LCS.store(use_lcs, Ordering::Relaxed);
}

pub(crate) fn use_lcs() -> bool {
    USE_LCS.load(Ordering::Relaxed)
}

/// The LCS array of `coloring` for a DBG of its SBWT, if the set operations use LCS arrays.
pub(crate) fn input_lcs<CSS: ColorSetStorage>(coloring: &CompactColexKmers<CSS>) -> Option<&LcsArray> {
    use_lcs().then(|| coloring.lcs()).filter(|lcs| !lcs.is_empty()) // Empty in an index built with --no-lcs
}

/// The DBG of the SBWT of `coloring`, an input of a set operation: from its LCS array if the set
/// operations use LCS arrays, and otherwise by column inversion, with 3-bit labels if
/// `optimize_peak_ram` (`--low-ram-mode`: less memory, more time, see Dbg::new_low_memory).
pub(crate) fn input_dbg<'a, CSS: ColorSetStorage>(coloring: &'a CompactColexKmers<CSS>, optimize_peak_ram: bool, n_threads: usize) -> Dbg<'a, SubsetMatrix> {
    match optimize_peak_ram {
        true => Dbg::new_low_memory(coloring.sbwt(), input_lcs(coloring), n_threads),
        false => Dbg::new(coloring.sbwt(), input_lcs(coloring), n_threads),
    }
}

pub(crate) fn empty_lcs() -> LcsArray {
    LcsArray::new(simple_sds_sbwt::int_vector::IntVector::default())
}
