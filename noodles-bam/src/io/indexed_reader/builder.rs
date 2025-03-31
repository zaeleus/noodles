use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, BinningIndex};

use super::IndexedReader;
use crate::bai;

/// An indexed BAM reader builder.
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
    /// use noodles_bam::{bai, io::indexed_reader::Builder};
    /// let index = bai::Index::default();
    /// let builder = Builder::default().set_index(index);
    /// ```
    pub fn set_index<I>(mut self, index: I) -> Self
    where
        I: BinningIndex + 'static,
    {
        self.index = Some(Box::new(index));
        self
    }

    /// Builds an indexed BAM reader from a path.
    ///
    /// If no index is set, this will attempt to read an associated index at `<src>.bai` or
    /// `<src>.csi`, in that order.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bam::io::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<bgzf::io::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let index = match self.index {
            Some(index) => index,
            None => read_associated_index(src)?,
        };

        let file = File::open(src)?;

        Ok(IndexedReader::new(file, index))
    }

    /// Builds an indexed BAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{bai, io::indexed_reader::Builder};
    /// let index = bai::Index::default();
    /// let data = [];
    /// let reader = Builder::default().set_index(index).build_from_reader(&data[..])?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<bgzf::io::Reader<R>>>
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
    let src = src.as_ref();

    match bai::fs::read(build_index_src(src, "bai")) {
        Ok(index) => Ok(Box::new(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let index = csi::fs::read(build_index_src(src, "csi"))?;
            Ok(Box::new(index))
        }
        Err(e) => Err(e),
    }
}

fn build_index_src<P, S>(src: P, ext: S) -> PathBuf
where
    P: AsRef<Path>,
    S: AsRef<OsStr>,
{
    push_ext(src.as_ref().into(), ext)
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
            push_ext(PathBuf::from("sample.bam"), "bai"),
            PathBuf::from("sample.bam.bai")
        );
    }
}
