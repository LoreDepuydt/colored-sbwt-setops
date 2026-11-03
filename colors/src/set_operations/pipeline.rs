//! The steps that all set operations on colored indexes share: interleaving the two SBWTs,
//! computing the SBWT of the result from the interleaving, and building the colors of the result
//! from the colors of the inputs.

use std::{cmp::max, sync::Arc};

use sbwt::{dbg::Dbg, LcsArray, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, util::u64_bitvec_to_simple_sds};
use crate::set_operations::{color_sets::{ColorCombination, ElementGeneratorFromInterleaving}, key_kmers::{mark_key_kmers_for, mark_key_kmers_from_boundaries, mark_structural_key_kmers, SubunitigBoundaries}, result_positions::InputToResult};
use crate::set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates};

/// The k-mers of the result of a set operation.
#[derive(Clone, Copy, Debug)]
pub(crate) enum SbwtOperation {
    /// The k-mers of both inputs. If keep_redundant_dummies is true, the dummy nodes that become
    /// redundant in the merge are kept: the result has the same k-mers and colors, but possibly
    /// more dummy nodes.
    Merge { keep_redundant_dummies: bool },
    /// The k-mers that are in both inputs.
    Intersect,
    /// The k-mers of the first input that are not in the second input.
    Difference,
}

impl SbwtOperation {
    /// How the result is called in the log
    fn result_name(&self) -> &'static str {
        match self {
            SbwtOperation::Merge { .. } => "merged",
            SbwtOperation::Intersect => "intersection",
            SbwtOperation::Difference => "difference",
        }
    }
}

type Bits = bitvec::vec::BitVec<u64, bitvec::order::Lsb0>;

/// Computes the result of `operation` on the k-mers of coloring1 and coloring2. The color set of
/// each k-mer of the result combines its color sets in the inputs as given by `colors`, and the
/// result has the colors `result_color_names`.
///
/// For the intersection and the difference, the SBWT of coloring1 needs select support.
#[allow(clippy::too_many_arguments)]
pub(crate) fn set_operation<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, operation: SbwtOperation, colors: ColorCombination, result_color_names: &[String], optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let interleaving = sbwt::MergeInterleaving::new(coloring1.sbwt(), coloring2.sbwt(), optimize_peak_ram, n_threads);
    assert_eq!(interleaving.s1.len(), interleaving.s2.len());

    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    // Temporarily destructure the colorings into parts in order to be able to put the SBWTs into
    // Arcs for the set operations of sbwt, which take Arcs by design (see the comment at sbwt::merge).
    let (sbwt1, lcs1, map1, sets1, color_names_1) = coloring1.into_parts();
    let (sbwt2, lcs2, map2, sets2, color_names_2) = coloring2.into_parts();
    // No k-mer of a difference is in coloring2, and its colors are those of coloring1, so the
    // difference needs only the SBWT of coloring2, for the set operation
    let is_difference = matches!(operation, SbwtOperation::Difference);
    let colors2 = if is_difference {
        drop((lcs2, map2, sets2, color_names_2));
        None
    } else {
        Some((lcs2, map2, sets2, color_names_2))
    };
    let sbwt1 = Arc::new(sbwt1);
    let sbwt2 = Arc::new(sbwt2);
    let interleaving = Arc::new(interleaving);

    // The result consists of the positions of an interleaving of the inputs that were not removed,
    // in order. Of that interleaving, only the positions of the inputs (s1 and s2) are needed
    // while dummy nodes have no colors, and the starts of its (k-1)-suffix groups, for the groups
    // of the result. The clones here clone just the Arcs.
    let (result_sbwt, s1, s2, is_leader, removed): (SbwtIndex<SubsetMatrix>, Bits, Bits, Bits, Option<Bits>) = match operation {
        SbwtOperation::Merge { keep_redundant_dummies } => {
            log::info!("Merging SBWTs");
            let (result_sbwt, removed) = if keep_redundant_dummies {
                (sbwt::merge_without_cleanup(sbwt1.clone(), sbwt2.clone(), interleaving.clone(), precalc_len, n_threads), None)
            } else {
                let (result_sbwt, removed) = sbwt::merge_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving.clone(), precalc_len, n_threads);
                (result_sbwt, Some(removed))
            };
            let sbwt::MergeInterleaving { s1, s2, is_leader, .. } = Arc::try_unwrap(interleaving).unwrap();
            (result_sbwt, s1, s2, is_leader, removed)
        },
        SbwtOperation::Intersect => {
            log::info!("Intersecting SBWTs");
            // The interleaving is passed on without keeping a reference, so that sbwt can free it
            // or take it over.
            let (result_sbwt, result_interleaving, removed) = sbwt::intersect_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving, precalc_len, optimize_peak_ram, n_threads)
                .unwrap_or_else(|_| panic!("The SBWT of the first index has no select support"));
            let sbwt::ThreeWayInterleaving { s1, s2, is_leader, .. } = result_interleaving;
            (result_sbwt, s1, s2, is_leader, Some(removed))
        },
        SbwtOperation::Difference => {
            log::info!("Subtracting SBWTs");
            // As for the intersection
            let (result_sbwt, result_interleaving, removed) = sbwt::difference_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving, precalc_len, optimize_peak_ram, n_threads)
                .unwrap_or_else(|_| panic!("The SBWT of the first index has no select support"));
            let sbwt::ThreeWayInterleaving { s1, s2, is_leader, .. } = result_interleaving;
            (result_sbwt, s1, s2, is_leader, Some(removed))
        },
    };

    // Put the coloring structs back together
    let sbwt1 = Arc::try_unwrap(sbwt1).unwrap();
    let sbwt2 = Arc::try_unwrap(sbwt2).unwrap();
    let coloring1 = CompactColexKmers::<CSS>::new(sbwt1, lcs1, map1, sets1, Some(&color_names_1));
    let coloring2 = match colors2 {
        Some((lcs2, map2, sets2, color_names_2)) => Some(CompactColexKmers::<CSS>::new(sbwt2, lcs2, map2, sets2, Some(&color_names_2))),
        None => {
            drop(sbwt2);
            None
        },
    };
    // Every k-mer of a merge is in one of the inputs, so it is visited when marking for the inputs
    // (and the check searches the result, which needs its LCS array)
    let check_visited = matches!(operation, SbwtOperation::Merge { .. }) && crate::set_operations::use_lcs();
    colored_result(&coloring1, None, coloring2.as_ref(), result_sbwt, &s1, &s2, is_leader, removed.as_ref(), operation.result_name(), colors, result_color_names, sample_distance, check_visited, optimize_peak_ram, n_threads)
}

/// Builds the colored index of the result of a set operation from its SBWT, which consists of the
/// positions of a merge plan of the inputs that were not removed, in order: `s1` and `s2` are the
/// positions of coloring1 and coloring2, and `removed` the removed positions (None if nothing was
/// removed). The color sets of the k-mers of the result combine their color sets in the inputs as
/// given by `colors`, and the result has the colors `result_color_names`. Without coloring2, its
/// k-mers are not used for the key k-mers, and `colors` must be ColorCombination::FromFirst.
/// `is_leader` are the starts of the (k-1)-suffix groups of the merge plan, from which the DBG of
/// the result takes its groups (see sbwt::result_k_minus_1_marks). `boundaries1` are the colored subunitigs of coloring1 if they are already known, so that its
/// unitigs are not walked again for the key k-mers; they are freed once the key k-mers are marked.
/// `check_visited` (see mark_key_kmers) needs that walk, so it can't be combined with `boundaries1`.
/// `result_name` is how the result is called in the log. With `optimize_peak_ram`, the DBGs of the
/// inputs that are built for the key k-mers take less memory and more time (see input_dbg).
#[allow(clippy::too_many_arguments)]
pub(crate) fn colored_result<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, boundaries1: Option<SubunitigBoundaries>, coloring2: Option<&CompactColexKmers<CSS>>, mut result_sbwt: SbwtIndex<SubsetMatrix>, s1: &Bits, s2: &Bits, is_leader: Bits, removed: Option<&Bits>, result_name: &str, colors: ColorCombination, result_color_names: &[String], sample_distance: usize, check_visited: bool, optimize_peak_ram: bool, n_threads: usize) -> CompactColexKmers<CSS> {
    // First, so that the merge plan's group starts are freed before anything else is built
    log::info!("Deriving the (k-1)-suffix groups of the {} SBWT from the merge plan", result_name);
    let k_minus_1_marks = sbwt::result_k_minus_1_marks(&is_leader, removed, n_threads);
    // The groups of the inputs whose unitigs mark_key_kmers walks, from the same merge plan, so
    // that their DBGs need no column inversion. Only for inputs without an LCS array to use (with
    // --no-lcs, or an index built without one): with it, the groups come from it at little cost,
    // and these would only take memory.
    let needs_marks = |coloring: &CompactColexKmers<CSS>| crate::set_operations::input_lcs(coloring).is_none();
    let marks1 = (boundaries1.is_none() && needs_marks(coloring1)).then(|| {
        log::info!("Deriving the (k-1)-suffix groups of the first input from the merge plan");
        input_k_minus_1_marks(&is_leader, s1, n_threads)
    });
    let marks2 = coloring2.filter(|c| needs_marks(c)).map(|_| {
        log::info!("Deriving the (k-1)-suffix groups of the second input from the merge plan");
        input_k_minus_1_marks(&is_leader, s2, n_threads)
    });
    let input_marks = (marks1, marks2);
    drop(is_leader);

    result_sbwt.build_select();

    let result_lcs = crate::set_operations::use_lcs().then(|| {
        log::info!("Building the LCS array for the {} SBWT", result_name);
        LcsArray::from_sbwt(&result_sbwt, n_threads, true)
    });

    log::info!("Initializing DBG for the {} SBWT", result_name);
    // Debug builds only (e.g. the tests): checks the groups from the merge plan against those from
    // the LCS array of the result, if it has one. Release builds skip this.
    if let (Some(lcs), true) = (result_lcs.as_ref(), cfg!(debug_assertions)) {
        assert_eq!(&k_minus_1_marks, Dbg::new(&result_sbwt, Some(lcs), n_threads).k_minus_1_marks(), "the groups of the {} from the merge plan", result_name);
    }
    let result_dbg = Dbg::new_with_k_minus_1_marks(&result_sbwt, k_minus_1_marks, n_threads);

    log::info!("=== Phase 1/3: marking key k-mers ===");
    let key_kmer_marks = mark_key_kmers(coloring1, boundaries1.as_ref(), coloring2, s1, s2, removed, input_marks, &result_sbwt, result_lcs.as_ref(), &result_dbg, sample_distance, check_visited, optimize_peak_ram, n_threads);
    // Only the key k-mers need these
    drop((result_dbg, boundaries1));
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / result_sbwt.n_kmers() as f64 * 100.0);

    let n_colors = u32::try_from(result_color_names.len()).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });
    let gen = || ElementGeneratorFromInterleaving {
        s1,
        s2,
        removed,
        coloring1,
        coloring2,
        result_key_kmer_marks: &key_kmer_marks,
        filter: None,
        colors,
        n_result_colors: n_colors as usize,
    };

    log::info!("=== PHASE 2/3: Building color set fingerprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    // TODO: do not clone key_kmer_marks
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id, _) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen(), key_kmer_marks.clone(), n_colors, n_threads, random_seed);

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let css = build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen(), n_threads);

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = crate::util::bitvec_to_simple_sds_bitvec(key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let colex_map = ColexToColorSetMap {
        sampling: key_kmer_marks,
        color_set_ids: key_kmer_idx_to_set_id,
    };

    let result_lcs = result_lcs.unwrap_or_else(crate::set_operations::empty_lcs);
    CompactColexKmers::<CSS>::new(result_sbwt, result_lcs, colex_map, css, Some(result_color_names))
}

/// The (k-1)-suffix groups of one input of a set operation, from the merge plan that the result
/// consists of: one bit per node of the input, set for the first node of each group, as
/// Dbg::k_minus_1_marks gives them. A group of the plan (from one `is_leader` bit to the next) is
/// the positions whose k-mers share their (k-1)-suffix, so the input's nodes in it form a group of
/// the input, which starts at the first position of the input in it (`in_input`: s1 or s2). New
/// dummy nodes of a dummy repair are in neither input, so they change nothing here.
///
/// Works on whole words of the plan, in parallel over chunks of words: a first pass counts the
/// nodes of the input in each chunk and finds whether its last group has a node of the input yet
/// (`started`), from which the starting node and `started` of every chunk follow; a second pass
/// marks. Within a word, each node of the input is a set bit of its word of `in_input`; it starts a
/// group of the input if no earlier node of the input is in its group of the plan, that is, if a
/// leader comes after the previous node of the input (in this word, or `started` before it).
fn input_k_minus_1_marks(is_leader: &Bits, in_input: &Bits, n_threads: usize) -> bitvec::vec::BitVec {
    use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
    assert_eq!(is_leader.len(), in_input.len());
    let n = in_input.len();
    let n_words = n.div_ceil(64);
    let (leaders, inputs) = (is_leader.as_raw_slice(), in_input.as_raw_slice());
    // Word w of a vector, without the bits after its end (a last word may have stale bits)
    let word = |v: &[u64], w: usize| if w + 1 == n_words && n % 64 != 0 { v[w] & ((1 << (n % 64)) - 1) } else { v[w] };
    // The bits after position p of a word (all of them if p is None)
    let after = |p: Option<u32>| match p { None => u64::MAX, Some(63) => 0, Some(p) => !((2u64 << p) - 1) };
    let chunk_words = n_words.div_ceil(n_threads.max(1) * 16).max(1);
    let chunks: Vec<std::ops::Range<usize>> = (0..n_words).step_by(chunk_words).map(|s| s..(s + chunk_words).min(n_words)).collect();
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    // Per chunk: its nodes of the input, whether it has a leader, and whether a node of the input
    // comes at or after its last leader (or anywhere in it, if it has none)
    let summaries: Vec<(usize, bool, bool)> = thread_pool.install(|| chunks.par_iter().map(|range| {
        let (mut nodes, mut has_leader, mut started) = (0, false, false);
        for w in range.clone() {
            let (l, s) = (word(leaders, w), word(inputs, w));
            nodes += s.count_ones() as usize;
            if l != 0 {
                has_leader = true;
                started = (s >> (63 - l.leading_zeros())) != 0;
            } else if s != 0 {
                started = true;
            }
        }
        (nodes, has_leader, started)
    }).collect());
    // The first node and `started` at the start of each chunk
    let mut starts = Vec::with_capacity(chunks.len());
    let (mut node, mut started) = (0, false);
    for &(nodes, has_leader, chunk_started) in &summaries {
        starts.push((node, started));
        node += nodes;
        started = if has_leader { chunk_started } else { started || chunk_started };
    }

    let marks = AtomicBitmap::new(node);
    thread_pool.install(|| chunks.par_iter().zip(starts).for_each(|(range, (mut node, mut started))| {
        for w in range.clone() {
            let (l, mut s) = (word(leaders, w), word(inputs, w));
            let mut prev = None; // the previous node of the input in this word
            while s != 0 {
                let b = s.trailing_zeros();
                let up_to_b = if b == 63 { u64::MAX } else { (2u64 << b) - 1 };
                if l & after(prev) & up_to_b != 0 {
                    started = false; // a group of the plan starts after the previous node of the input
                }
                if !started {
                    marks.set(node, true);
                    started = true;
                }
                node += 1;
                prev = Some(b);
                s &= s - 1;
            }
            if l & after(prev) != 0 {
                started = false;
            }
        }
    }));
    marks.into_bitvec()
}

/// Marks the key k-mers of the result: the structural key k-mers of the result, and the k-mers of
/// the result where a colored subunitig of an input ends (see mark_key_kmers_for). In an
/// intersection, a color set changes between consecutive k-mers x and y of a unitig of the result
/// only if the color set of x or y changes in one of the inputs, which has the edge from x to y
/// too. Then x is the last k-mer of a colored subunitig of that input, so it is marked. In a
/// difference, the color sets are those of coloring1, which has the edge from x to y too, so
/// marking for coloring1 suffices, and coloring2 is None.
///
/// With `boundaries1`, the key k-mers for coloring1 are marked from its known colored subunitigs
/// instead of a walk over its unitigs. `input_marks` are the (k-1)-suffix groups of the walked
/// inputs (coloring1 if it is walked, and coloring2), if known (see input_k_minus_1_marks); the
/// DBGs of the others compute them.
///
/// With check_visited, checks in debug builds that every k-mer of the result is visited while
/// marking for the inputs (see mark_key_kmers_for). Skipped in release builds: it costs a full
/// extra streaming-index pass over the result and is not needed for marking key k-mers, only for
/// catching bugs during development.
#[allow(clippy::too_many_arguments)]
fn mark_key_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, boundaries1: Option<&SubunitigBoundaries>, coloring2: Option<&CompactColexKmers<CSS>>, s1: &Bits, s2: &Bits, removed: Option<&Bits>, input_marks: (Option<bitvec::vec::BitVec>, Option<bitvec::vec::BitVec>), result_sbwt: &SbwtIndex<SubsetMatrix>, result_lcs: Option<&LcsArray>, result_dbg: &Dbg<'_, SubsetMatrix>, sample_distance: usize, check_visited: bool, optimize_peak_ram: bool, n_threads: usize) -> bitvec::vec::BitVec {
    let k = result_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert!(coloring2.is_none_or(|coloring2| coloring2.get_k() == k));

    let key_kmer_marks = AtomicBitmap::new(result_sbwt.n_sets());
    mark_structural_key_kmers(result_sbwt, result_dbg, &key_kmer_marks, sample_distance, n_threads);

    let visited_marks = (check_visited && cfg!(debug_assertions)).then(|| AtomicBitmap::new(result_sbwt.n_sets()));
    // TODO: This is a second copy of the removed positions (n bits), next to the original that
    // colors 2 and 3 use. Converting them to simple-sds once, with rank support, and sharing that
    // with colors 2 and 3 would save the copy.
    let removed = removed.map(|r| {
        let mut r = u64_bitvec_to_simple_sds(r);
        r.enable_rank();
        r
    });
    if let Some(boundaries1) = boundaries1 {
        let to_result = InputToResult::new(s1, removed.as_ref());
        mark_key_kmers_from_boundaries(boundaries1, |colex| to_result.result_colex(colex), result_dbg, &key_kmer_marks, n_threads);
    }
    let (marks1, marks2) = input_marks;
    let walked_coloring1 = boundaries1.is_none().then_some((coloring1, s1, marks1));
    for (coloring, in_input, marks) in walked_coloring1.into_iter().chain(coloring2.map(|coloring2| (coloring2, s2, marks2))) {
        let to_result = InputToResult::new(in_input, removed.as_ref());
        mark_key_kmers_for(coloring, |colex| to_result.result_colex(colex), marks, result_sbwt, result_lcs, result_dbg, &key_kmer_marks, visited_marks.as_ref(), optimize_peak_ram, n_threads);
    }
    drop(removed);

    if let Some(visited_marks) = visited_marks {
        assert_eq!(visited_marks.into_bitvec().count_ones(), result_sbwt.n_kmers());
    }

    key_kmer_marks.into_bitvec()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use sbwt::dbg::Dbg;

    use crate::set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING};

    /// The groups of each input from the merge plans of a merge, an intersection, a difference and a
    /// filter (with dummy repairs) against the groups its own DBG computes from its SBWT.
    #[test]
    fn input_k_minus_1_marks_from_plans() {
        for (k, seed) in [(3, 1), (5, 2), (7, 3), (11, 4), (31, 5)] {
            // Overlapping inputs: the second has a part of the first, so the plans interleave
            let a = gen_random_dna_string(3000, seed);
            let b: Vec<u8> = a[1000..2000].iter().copied().chain(gen_random_dna_string(1500, seed + 100)).collect();
            let index = |seqs: Vec<(&str, Vec<u8>)>| {
                let mut names = Vec::<String>::new();
                let seqs = assign_color_ids(&mut names, &seqs);
                build_named_coloring(k, &names, &seqs, NO_SAMPLING, 2)
            };
            let (c1, c2) = (index(vec![("a", a.clone())]), index(vec![("b", b.clone())]));
            let (sbwt1, sbwt2) = (c1.sbwt(), c2.sbwt());
            let groups1 = Dbg::new(sbwt1, None, 2).k_minus_1_marks().clone();
            let groups2 = Dbg::new(sbwt2, None, 2).k_minus_1_marks().clone();
            let check = |what: &str, is_leader: &super::Bits, s1: &super::Bits, s2: Option<&super::Bits>| {
                for n_threads in [1, 3, 16] {
                    assert_eq!(super::input_k_minus_1_marks(is_leader, s1, n_threads), groups1, "k={k}: {what}, first input, {n_threads} threads");
                    if let Some(s2) = s2 {
                        assert_eq!(super::input_k_minus_1_marks(is_leader, s2, n_threads), groups2, "k={k}: {what}, second input, {n_threads} threads");
                    }
                }
            };

            let merge = sbwt::MergeInterleaving::new(sbwt1, sbwt2, false, 2);
            check("merge", &merge.is_leader, &merge.s1, Some(&merge.s2));
            let (arc1, arc2) = (Arc::new(sbwt1.clone()), Arc::new(sbwt2.clone()));
            let (_, inter, _) = sbwt::intersect_with_removed_positions(arc1.clone(), arc2.clone(), Arc::new(sbwt::MergeInterleaving::new(sbwt1, sbwt2, false, 2)), 0, false, 2).unwrap();
            check("intersection", &inter.is_leader, &inter.s1, Some(&inter.s2));
            let (_, diff, _) = sbwt::difference_with_removed_positions(arc1.clone(), arc2, Arc::new(sbwt::MergeInterleaving::new(sbwt1, sbwt2, false, 2)), 0, false, 2).unwrap();
            check("difference", &diff.is_leader, &diff.s1, Some(&diff.s2));
            // A filter keeping every third k-mer: many dummy repairs (s2 is the kept nodes, not an input)
            let keep: super::Bits = (0..sbwt1.n_sets()).map(|i| i % 3 == 0).collect();
            let (_, filtered, _) = sbwt::filter_with_removed_positions(arc1, keep, None, 0, false, 2).unwrap();
            check("filter", &filtered.is_leader, &filtered.s1, None);
        }
    }

    /// The marks bit by bit, as the definition says: the reference for the word-wise version.
    fn input_k_minus_1_marks_by_bit(is_leader: &super::Bits, in_input: &super::Bits) -> bitvec::vec::BitVec {
        let mut marks = bitvec::bitvec![0; in_input.count_ones()];
        let (mut node, mut started) = (0, false);
        for pos in 0..in_input.len() {
            if is_leader[pos] {
                started = false;
            }
            if in_input[pos] {
                if !started {
                    marks.set(node, true);
                    started = true;
                }
                node += 1;
            }
        }
        marks
    }

    /// The word-wise version against the bit-by-bit one on random plans: lengths around words and
    /// chunks, densities from sparse to full, stale bits after the end, and thread counts that put
    /// chunk borders everywhere.
    #[test]
    fn input_k_minus_1_marks_word_wise() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(11);
        for len in [0, 1, 2, 63, 64, 65, 127, 128, 129, 1000, 4096, 70001] {
            for (leader_pct, input_pct) in [(0, 50), (1, 50), (10, 10), (50, 90), (100, 100), (5, 0), (30, 3)] {
                let mut bits = |pct: u32| -> super::Bits { (0..len).map(|_| rng.next_u32() % 100 < pct).collect() };
                let (mut is_leader, mut in_input) = (bits(leader_pct), bits(input_pct));
                if len > 0 { is_leader.set(0, true); } // as in a plan: position 0 starts a group
                // Stale bits after the end, as after a truncation
                for v in [&mut is_leader, &mut in_input] {
                    v.extend(std::iter::repeat_n(true, 70));
                    v.truncate(len);
                }
                let expected = input_k_minus_1_marks_by_bit(&is_leader, &in_input);
                for n_threads in [1, 2, 7, 64] {
                    assert_eq!(super::input_k_minus_1_marks(&is_leader, &in_input, n_threads), expected, "len {len}, leaders {leader_pct} %, input {input_pct} %, {n_threads} threads");
                }
            }
        }
    }
}
