//! Finding the k-mers of a sequence in an SBWT, with or without its LCS array.
//!
//! With an LCS array, the k-mers are streamed with matching statistics. An SBWT built without one
//! (`build --no-lcs`) finds them by following the edges of the de Bruijn graph instead, using its
//! (k-1)-suffix groups, which are computed from the SBWT alone (column inversion).

use std::ops::Range;

use bitvec::vec::BitVec;
use sbwt::{ExtendRight, LcsArray, MatchingStatisticsIterator, SbwtIndex, StreamingIndex, SubsetMatrix};

/// Returns the LCS array, or None if it is empty (an SBWT built without one).
pub fn lcs_if_present(lcs: &LcsArray) -> Option<&LcsArray> {
    (!lcs.is_empty()).then_some(lcs)
}

/// Finds k-mers in an SBWT, through its LCS array or through its (k-1)-suffix groups.
#[derive(Clone, Copy)]
pub enum KmerLocator<'a> {
    Lcs { sbwt: &'a SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray },
    SuffixGroups { sbwt: &'a SbwtIndex<SubsetMatrix>, k_minus_1_marks: &'a BitVec },
}

impl<'a> KmerLocator<'a> {
    /// Uses the LCS array if it is not empty. Otherwise `k_minus_1_marks` must be given.
    pub fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, k_minus_1_marks: Option<&'a BitVec>) -> Self {
        match (lcs_if_present(lcs), k_minus_1_marks) {
            (Some(lcs), _) => Self::Lcs { sbwt, lcs },
            (None, Some(k_minus_1_marks)) => Self::SuffixGroups { sbwt, k_minus_1_marks },
            (None, None) => panic!("An SBWT without an LCS array needs its (k-1)-suffix groups to find k-mers"),
        }
    }

    pub fn sbwt(&self) -> &'a SbwtIndex<SubsetMatrix> {
        match *self {
            Self::Lcs { sbwt, .. } | Self::SuffixGroups { sbwt, .. } => sbwt,
        }
    }

    pub fn k(&self) -> usize {
        self.sbwt().k()
    }

    /// Like [StreamingIndex::matching_statistics_iter]: one (match length, colex range) pair per
    /// position of `query`. Without an LCS array, only k-mers are matched: the length is k if the
    /// k-mer ending at the position is in the SBWT, and 0 otherwise.
    pub fn matching_statistics_iter<'q>(&self, query: &'q [u8]) -> KmerMatches<'a, 'q> {
        match *self {
            Self::Lcs { sbwt, lcs } => KmerMatches::Lcs(StreamingIndex::new(sbwt, lcs).matching_statistics_iter(query)),
            Self::SuffixGroups { sbwt, k_minus_1_marks } => KmerMatches::SuffixGroups(SuffixGroupMatches {
                sbwt, k_minus_1_marks, query, pos: 0, prev: None, last_invalid: None,
            }),
        }
    }
}

pub enum KmerMatches<'a, 'q> {
    Lcs(MatchingStatisticsIterator<'a, 'q, SbwtIndex<SubsetMatrix>, LcsArray>),
    SuffixGroups(SuffixGroupMatches<'a, 'q>),
}

impl Iterator for KmerMatches<'_, '_> {
    type Item = (usize, Range<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Lcs(iter) => iter.next(),
            Self::SuffixGroups(iter) => iter.next(),
        }
    }
}

/// Finds the k-mers of a query through the (k-1)-suffix groups of the SBWT: the k-mer after a
/// found k-mer is its out-neighbor, one SBWT step from its suffix group. Other k-mers are searched
/// from scratch.
pub struct SuffixGroupMatches<'a, 'q> {
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    k_minus_1_marks: &'a BitVec,
    query: &'q [u8],
    pos: usize,
    prev: Option<usize>, // Colex rank of the k-mer ending at pos-1, if found
    last_invalid: Option<usize>, // Position of the last non-ACGT character so far
}

impl Iterator for SuffixGroupMatches<'_, '_> {
    type Item = (usize, Range<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        let end = self.pos;
        if end >= self.query.len() {
            return None;
        }
        self.pos += 1;
        let k = self.sbwt.k();
        if !is_acgt(self.query[end]) {
            self.last_invalid = Some(end);
        }
        if end + 1 < k {
            return Some((0, 0..0));
        }
        let start = end + 1 - k;
        let cur = if self.last_invalid.is_some_and(|p| p >= start) {
            None
        } else if let Some(v) = self.prev {
            follow_outedge(self.sbwt, self.k_minus_1_marks, v, self.query[end])
        } else {
            self.sbwt.search(&self.query[start..=end]).map(|range| range.start)
        };
        self.prev = cur;
        Some(match cur {
            Some(v) => (k, v..v + 1),
            None => (0, 0..0),
        })
    }
}

/// The colex rank of the out-neighbor of node `v` with label `c`, if it exists.
fn follow_outedge(sbwt: &SbwtIndex<SubsetMatrix>, k_minus_1_marks: &BitVec, v: usize, c: u8) -> Option<usize> {
    let mut group_start = v;
    while !k_minus_1_marks[group_start] {
        group_start -= 1; // Terminates because the first node starts a group
    }
    let mut group_end = v + 1;
    while group_end < k_minus_1_marks.len() && !k_minus_1_marks[group_end] {
        group_end += 1;
    }
    let range = sbwt.extend_right(group_start..group_end, c);
    (!range.is_empty()).then_some(range.start)
}

pub(crate) fn is_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use sbwt::{dbg::Dbg, BitPackedKmerSortingMem};

    // Without an LCS array, the k-mers found must be exactly those found by searching every k-mer,
    // including around mismatches, Ns and lowercase characters.
    #[test]
    fn suffix_groups_find_the_kmers() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let random_seq = |rng: &mut rand::rngs::StdRng, len: usize| -> Vec<u8> { (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect() };
        for k in [3, 5, 15, 31] {
            let seqs: Vec<Vec<u8>> = (0..20).map(|_| random_seq(&mut rng, 200)).collect();
            let slices: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
            let (mut sbwt, lcs) = BitPackedKmerSortingMem::new_from_slices(&slices, k).add_rev_comp(true).build_lcs(true).run();
            sbwt.build_select();
            let marks = Dbg::new(&sbwt, None, 2).k_minus_1_marks().clone(); // Without the LCS array
            assert_eq!(&marks, Dbg::new(&sbwt, lcs.as_ref(), 2).k_minus_1_marks());
            let empty = crate::set_operations::empty_lcs();
            let locator = KmerLocator::new(&sbwt, &empty, Some(&marks));

            let mut queries = seqs.clone();
            queries.extend((0..20).map(|_| random_seq(&mut rng, 100)));
            for seq in &seqs[..5] {
                let mut q = seq.clone();
                q[50] = b'N';
                q[120] = if q[120] == b'A' { b'C' } else { b'A' };
                q[150] = q[150].to_ascii_lowercase();
                queries.push(q);
            }
            for q in &queries {
                let found: Vec<Option<usize>> = locator.matching_statistics_iter(q).skip(k - 1)
                    .map(|(len, range)| (len == k).then_some(range.start)).collect();
                let searched: Vec<Option<usize>> = q.windows(k)
                    .map(|w| if w.iter().all(|&c| is_acgt(c)) { sbwt.search(w).map(|r| r.start) } else { None }).collect();
                assert_eq!(found, searched, "k = {}", k);
            }
        }
    }
}
