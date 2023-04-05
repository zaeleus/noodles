use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_csi as csi;

use super::IndexedReader;
use crate::bai;

/// An indexed BAM reader builder.
#[derive(Default)]
pub struct Builder {
    index: Option<csi::Index>,
}

impl Builder {
    /// Sets an index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::indexed_reader::Builder;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let builder = Builder::default().set_index(index);
    /// ```
    pub fn set_index(mut self, index: csi::Index) -> Self {
        self.index = Some(index);
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
    /// use noodles_bam::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<IndexedReader<bgzf::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.index.is_none() {
            self.index = read_associated_index(src).map(Some)?;
        }

        let file = File::open(src)?;
        self.build_from_reader(file)
    }

    /// Builds an indexed BAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::indexed_reader::Builder;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let data = [];
    /// let builder = Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(&data[..])?;
    /// # Ok::<_, std::io::Error>(())
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

fn read_associated_index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();

    match bai::read(build_index_src(src, "bai")) {
        Ok(index) => Ok(index),
        Err(e) if e.kind() == io::ErrorKind::NotFound => csi::read(build_index_src(src, "csi")),
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
