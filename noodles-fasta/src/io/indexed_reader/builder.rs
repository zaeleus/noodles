use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;

use super::IndexedReader;
use crate::fai;

/// An indexed FASTA reader builder.
#[derive(Default)]
pub struct Builder {
    index: Option<fai::Index>,
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
            Some("gz" | "bgz") => bgzf::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(crate::io::BufReader::Bgzf)?,
            _ => File::open(src)
                .map(BufReader::new)
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
}
