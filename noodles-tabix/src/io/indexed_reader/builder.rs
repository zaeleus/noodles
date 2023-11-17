use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io,
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, io::IndexedReader};

/// An indexed reader builder.
#[derive(Default)]
pub struct Builder {
    index: Option<csi::Index>,
}

impl Builder {
    /// Sets an index.
    pub fn set_index(mut self, index: csi::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed reader from a path.
    pub fn build_from_path<P>(
        self,
        src: P,
    ) -> io::Result<IndexedReader<bgzf::Reader<File>, csi::Index>>
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
}

fn read_associated_index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    crate::read(build_index_src(src, "tbi"))
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
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_push_ext() {
        assert_eq!(
            push_ext(PathBuf::from("annotations.gff.gz"), "tbi"),
            PathBuf::from("annotations.gff.gz.tbi")
        );
    }
}
