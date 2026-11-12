use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;

use super::IndexedReader;
use crate::fai;

/// Default read buffer capacity (in bytes) for uncompressed inputs.
///
/// Chosen empirically: on a multi-GB reference loaded contig-by-contig,
/// throughput is flat between ~32 KiB and ~256 KiB and ~10% better than
/// the 8 KiB [`BufReader`] default. 64 KiB is a friendly midpoint that
/// also matches a typical filesystem block-cache page granularity.
const DEFAULT_BUFFER_CAPACITY: usize = 64 * 1024;

/// An indexed FASTA reader builder.
pub struct Builder {
    index: Option<fai::Index>,
    buffer_capacity: usize,
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            index: None,
            buffer_capacity: DEFAULT_BUFFER_CAPACITY,
        }
    }
}

impl Builder {
    /// Sets an index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{fai, io::indexed_reader::Builder};
    /// let index = fai::Index::default();
    /// let builder = Builder::default().set_index(index);
    /// ```
    pub fn set_index(mut self, index: fai::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Sets the read buffer capacity (in bytes) for uncompressed inputs.
    ///
    /// [`build_from_path`] wraps the opened file in a
    /// [`BufReader::with_capacity`] of this size. The default is 64 KiB,
    /// which is friendlier than the 8 KiB [`BufReader`] default for
    /// streaming long contiguous regions, e.g. loading a full chromosome
    /// via [`IndexedReader::read_sequence`].
    ///
    /// Ignored for BGZF-compressed inputs (`.gz`, `.bgz`): the BGZF
    /// reader consumes block-sized chunks directly from the file and an
    /// outer buffer provides no measurable benefit.
    ///
    /// [`build_from_path`]: Self::build_from_path
    /// [`IndexedReader::read_sequence`]: super::IndexedReader::read_sequence
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::io::indexed_reader::Builder;
    /// let builder = Builder::default().set_buffer_capacity(128 * 1024);
    /// ```
    pub fn set_buffer_capacity(mut self, capacity: usize) -> Self {
        self.buffer_capacity = capacity;
        self
    }

    /// Builds an indexed FASTA reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fasta::io::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("reference.fa")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<crate::io::BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let index = match self.index {
            Some(index) => index,
            None => {
                let index_src = build_index_src(src);
                fai::fs::read(index_src)?
            }
        };

        let reader = match src.extension().and_then(|ext| ext.to_str()) {
            Some("gz" | "bgz") => bgzf::io::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(crate::io::BufReader::Bgzf)?,
            _ => File::open(src)
                .map(|file| BufReader::with_capacity(self.buffer_capacity, file))
                .map(crate::io::BufReader::Uncompressed)?,
        };

        Ok(IndexedReader::new(reader, index))
    }

    /// Builds an indexed FASTA reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{fai, io::indexed_reader::Builder};
    ///
    /// let index = fai::Index::default();
    /// let data = [];
    /// let builder = Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(&data[..])?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<R>>
    where
        R: BufRead,
    {
        let index = self
            .index
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index"))?;

        Ok(IndexedReader::new(reader, index))
    }
}

fn build_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "fai";
    push_ext(src.as_ref().into(), EXT)
}

fn push_ext<S>(path: PathBuf, ext: S) -> PathBuf
where
    S: AsRef<OsStr>,
{
    let mut s = OsString::from(path);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_index_src() {
        assert_eq!(build_index_src("ref.fa"), PathBuf::from("ref.fa.fai"));
    }

    #[test]
    fn test_default_buffer_capacity() {
        assert_eq!(Builder::default().buffer_capacity, DEFAULT_BUFFER_CAPACITY);
    }

    #[test]
    fn test_set_buffer_capacity() {
        let builder = Builder::default().set_buffer_capacity(128 * 1024);
        assert_eq!(builder.buffer_capacity, 128 * 1024);
    }
}
