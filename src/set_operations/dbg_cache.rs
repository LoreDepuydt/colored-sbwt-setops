//! For development only: the marks of the DBG of the input of a subset, stored on disk, so that
//! later subsets of the same index skip computing them (`themisto2 subset --dbg-cache`). They are
//! its (k-1)-suffix groups (from the column inversion without an LCS) and its dummy nodes, one bit
//! per node each, and depend only on the SBWT of the index, not on the kept colors.

use std::{fs::File, io::{BufReader, BufWriter, ErrorKind, Read, Write}, path::{Path, PathBuf}, time::UNIX_EPOCH};

use sbwt::{dbg::Dbg, SbwtIndex, SubsetMatrix};

type Marks = bitvec::vec::BitVec;

// The marks are stored as their words, which are usize
const _: () = assert!(usize::BITS == 64);

const MAGIC: &[u8; 16] = b"THEMISTO2DBGMRK1";

/// The stored DBG marks of one index: a file in the cache directory named after the index file.
/// The file starts with what the marks belong to (the path, size and modification time of the
/// index file, and the number of nodes and k of its SBWT), and is used only if all of that matches.
pub(crate) struct DbgCache {
    file: PathBuf,
    index_key: String,
}

impl DbgCache {
    pub(crate) fn new(dir: &Path, index: &Path) -> Self {
        let path = std::fs::canonicalize(index).unwrap_or_else(|_| index.to_path_buf());
        let meta = std::fs::metadata(&path).unwrap_or_else(|e| panic!("Could not read the metadata of {}: {}", path.display(), e));
        let mtime = meta.modified().ok().and_then(|t| t.duration_since(UNIX_EPOCH).ok()).map_or(0, |d| d.as_nanos());
        let name = path.file_name().map_or("index".into(), |n| n.to_string_lossy().into_owned());
        Self { file: dir.join(format!("{name}.dbg-marks")), index_key: format!("{}\t{}\t{}", path.display(), meta.len(), mtime) }
    }

    /// The DBG of `sbwt` (the SBWT of the index), with its marks from the cache if it has them for
    /// this index, and otherwise computed by `compute` and then stored. Failing to store them is
    /// logged and does not stop the subset.
    pub(crate) fn dbg<'a>(&self, sbwt: &'a SbwtIndex<SubsetMatrix>, compute: impl FnOnce() -> Dbg<'a, SubsetMatrix>) -> Dbg<'a, SubsetMatrix> {
        match self.load(sbwt) {
            Ok(Some((k_minus_1_marks, dummy_marks))) => {
                log::info!("Took the DBG marks from {}: skipping their computation", self.file.display());
                return Dbg::from_marks(sbwt, k_minus_1_marks, dummy_marks);
            },
            Ok(None) => log::info!("No DBG marks for this index in {}: computing them", self.file.display()),
            Err(e) => log::warn!("Could not read the DBG marks from {}: {}; computing them", self.file.display(), e),
        }
        let dbg = compute();
        log::info!("Storing the DBG marks in {}", self.file.display());
        if let Err(e) = self.store(sbwt, &dbg) {
            log::warn!("Could not store the DBG marks in {}: {}", self.file.display(), e);
        }
        dbg
    }

    /// What the marks belong to, as stored at the start of the file
    fn header(&self, sbwt: &SbwtIndex<SubsetMatrix>) -> String {
        format!("{}\t{}\t{}\n", self.index_key, sbwt.n_sets(), sbwt.k())
    }

    /// The (k-1)-suffix group marks and dummy marks, or None if the cache has none for this index.
    fn load(&self, sbwt: &SbwtIndex<SubsetMatrix>) -> std::io::Result<Option<(Marks, Marks)>> {
        let file = match File::open(&self.file) {
            Err(e) if e.kind() == ErrorKind::NotFound => return Ok(None),
            file => file?,
        };
        let mut input = BufReader::with_capacity(1 << 23, file);
        let mut magic = [0u8; 16];
        input.read_exact(&mut magic)?;
        let mut header_len = [0u8; 8];
        input.read_exact(&mut header_len)?;
        let mut header = vec![0u8; u64::from_le_bytes(header_len) as usize];
        input.read_exact(&mut header)?;
        if &magic != MAGIC || header != self.header(sbwt).as_bytes() {
            log::info!("The DBG marks in {} are for another index: {}", self.file.display(), String::from_utf8_lossy(&header).trim_end());
            return Ok(None);
        }
        let k_minus_1_marks = read_marks(&mut input, sbwt.n_sets())?;
        let dummy_marks = read_marks(&mut input, sbwt.n_sets())?;
        Ok(Some((k_minus_1_marks, dummy_marks)))
    }

    /// Writes a temporary file and renames it, so that an interrupted write leaves no cache.
    fn store(&self, sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<'_, SubsetMatrix>) -> std::io::Result<()> {
        let tmp = self.file.with_extension("dbg-marks.tmp");
        let mut out = BufWriter::with_capacity(1 << 23, File::create(&tmp)?);
        let header = self.header(sbwt);
        out.write_all(MAGIC)?;
        out.write_all(&(header.len() as u64).to_le_bytes())?;
        out.write_all(header.as_bytes())?;
        for marks in [dbg.k_minus_1_marks(), dbg.dummy_marks()] {
            for &word in &marks.as_raw_slice()[..marks.len().div_ceil(64)] {
                out.write_all(&(word as u64).to_le_bytes())?;
            }
        }
        // Flushed first, so that a failed write is reported rather than renamed
        out.into_inner().map_err(|e| e.into_error())?;
        std::fs::rename(&tmp, &self.file)
    }
}

/// `n` bits stored as their little-endian 64-bit words.
fn read_marks(input: &mut impl Read, n: usize) -> std::io::Result<Marks> {
    let mut words = vec![0usize; n.div_ceil(64)];
    let mut buf = vec![0u8; 1 << 23];
    for chunk in words.chunks_mut(buf.len() / 8) {
        let bytes = &mut buf[..chunk.len() * 8];
        input.read_exact(bytes)?;
        for (word, b) in chunk.iter_mut().zip(bytes.chunks_exact(8)) {
            *word = u64::from_le_bytes(b.try_into().unwrap()) as usize;
        }
    }
    let mut marks = Marks::from_vec(words);
    marks.truncate(n);
    Ok(marks)
}

#[cfg(test)]
mod tests {
    use sbwt::dbg::Dbg;

    use crate::set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING};

    /// Computes the DBG marks of a small index through a cache, takes them from the cache, and
    /// checks that they are the same as computed directly; then checks that a changed index file
    /// does not take them. Uses a directory under $TMPDIR.
    #[test]
    fn store_and_take() {
        let dir = std::env::temp_dir().join(format!("themisto2_dbg_cache_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let index = dir.join("index.thm2");
        std::fs::write(&index, b"stands for the index file").unwrap();

        let mut names = Vec::<String>::new();
        let colors = [("a", gen_random_dna_string(3000, 1)), ("b", gen_random_dna_string(2000, 2))];
        let seqs = assign_color_ids(&mut names, &colors);
        let coloring = build_named_coloring(7, &names, &seqs, NO_SAMPLING, 3);
        let sbwt = coloring.sbwt();
        let direct = Dbg::new(sbwt, None, 3);

        let cache = super::DbgCache::new(&dir, &index);
        let file = dir.join("index.thm2.dbg-marks");
        assert!(cache.load(sbwt).unwrap().is_none());
        let computed = cache.dbg(sbwt, || Dbg::new(sbwt, None, 3));
        assert!(file.exists());
        let (k_minus_1_marks, dummy_marks) = cache.load(sbwt).unwrap().expect("the marks are stored");
        for dbg in [&computed, &cache.dbg(sbwt, || Dbg::new(sbwt, None, 3))] {
            assert_eq!(dbg.k_minus_1_marks(), direct.k_minus_1_marks());
            assert_eq!(dbg.dummy_marks(), direct.dummy_marks());
        }
        assert_eq!(&k_minus_1_marks, direct.k_minus_1_marks());
        assert_eq!(&dummy_marks, direct.dummy_marks());

        // Another size of the index file: the stored marks are not for it
        std::fs::write(&index, b"stands for another index file").unwrap();
        assert!(super::DbgCache::new(&dir, &index).load(sbwt).unwrap().is_none());

        std::fs::remove_dir_all(&dir).unwrap();
    }
}
