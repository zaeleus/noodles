use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io,
    path::{Path, PathBuf},
};

use super::IndexedReader;
use crate::{fai, reader::Source};

/// An indexed reader adapter builder.
///
/// This is a convenience builder for creating an indexed reader adapter from paths on a
/// filesystem.
///
/// By default, it opens a [`Reader`] for a source path (`src`) and reads its associated index at
/// `<src>.fai`. The index location can be overridden by calling [`set_index_src`].
#[derive(Default)]
pub struct Builder {
    index_src: Option<PathBuf>,
}

impl Builder {
    /// Sets the index source path.
    ///
    /// When set, this path is used instead of inferring one from the given source path.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::repository::adapters::IndexedReader;
    /// let builder = IndexedReader::builder().set_index_src("reference.fa.fai");
    /// ```
    pub fn set_index_src<P>(mut self, index_src: P) -> Self
    where
        P: Into<PathBuf>,
    {
        self.index_src = Some(index_src.into());
        self
    }

    /// Creates an indexed reader adapter from the given path.
    ///
    /// By default, `<src>.fai` is used as the path to the associated index. This can be overridden
    /// by calling [`set_index_src`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_fasta::repository::adapters::IndexedReader;
    /// let adapter = IndexedReader::builder().open("reference.fa")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn open<P>(self, src: P) -> io::Result<IndexedReader<Source<File>>>
    where
        P: AsRef<Path>,
    {
        use crate::reader::Builder;

        let index_src = self
            .index_src
            .unwrap_or_else(|| push_ext(src.as_ref().to_path_buf(), "fai"));
        let index = fai::read(index_src)?;

        let reader = Builder::default()
            .set_fasta_index(index)
            .build_from_path(src)?;

        Ok(IndexedReader::new(reader, Default::default()))
    }
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
    fn test_push_ext() {
        assert_eq!(
            push_ext(PathBuf::from("reference.fa"), "fai"),
            PathBuf::from("reference.fa.fai")
        );
    }
}
