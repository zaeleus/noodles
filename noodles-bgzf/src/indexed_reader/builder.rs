use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use super::IndexedReader;
use crate::{gzi, reader};

/// An indexed BGZF reader builder.
#[derive(Default)]
pub struct Builder {
    reader_builder: reader::Builder,
    index: Option<gzi::Index>,
}

impl Builder {
    /// Sets a GZ index.
    pub fn set_index(mut self, index: gzi::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed BGZF reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let index = match self.index {
            Some(index) => index,
            None => {
                let index_src = build_index_src(&src);
                gzi::fs::read(index_src)?
            }
        };

        let inner = self.reader_builder.build_from_path(src)?;

        Ok(IndexedReader { inner, index })
    }

    /// Builds a indexed BGZF reader from a reader.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<R>>
    where
        R: Read,
    {
        let index = self
            .index
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index"))?;

        let inner = self.reader_builder.build_from_reader(reader);

        Ok(IndexedReader { inner, index })
    }
}

fn build_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "gzi";
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
            build_index_src("noodles.gz"),
            PathBuf::from("noodles.gz.gzi")
        );
    }
}
