use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, BinningIndex};

use super::IndexedReader;

/// An indexed BCF reader.
#[derive(Default)]
pub struct Builder {
    index: Option<Box<dyn BinningIndex>>,
}

impl Builder {
    /// Sets an index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::io::indexed_reader::Builder;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let builder = Builder::default().set_index(index);
    /// ```
    pub fn set_index<I>(mut self, index: I) -> Self
    where
        I: BinningIndex + 'static,
    {
        self.index = Some(Box::new(index));
        self
    }

    /// Builds an indexed BCF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::io::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bcf")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<bgzf::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let file = File::open(src)?;

        let index = match self.index {
            Some(index) => index,
            None => read_associated_index(src)?,
        };

        Ok(IndexedReader::new(file, index))
    }

    /// Builds an indexed BCF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::io::indexed_reader::Builder;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let reader = Builder::default().set_index(index).build_from_reader(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<bgzf::Reader<R>>>
    where
        R: Read,
    {
        let index = self
            .index
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index"))?;

        Ok(IndexedReader::new(reader, index))
    }
}

fn read_associated_index<P>(src: P) -> io::Result<Box<dyn BinningIndex>>
where
    P: AsRef<Path>,
{
    let index = csi::fs::read(build_index_src(src))?;
    Ok(Box::new(index))
}

fn build_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "csi";
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
        assert_eq!(
            build_index_src("sample.bcf"),
            PathBuf::from("sample.bcf.csi")
        );
    }
}
