use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use noodles_csi::{self as csi, BinningIndex};

use super::IndexedReader;

/// An indexed SAM reader builder.
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
    /// use noodles_csi as csi;
    /// use noodles_sam::indexed_reader::Builder;
    ///
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

    /// Builds an indexed SAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_sam::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.sam.gz")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.index.is_none() {
            let index_src = build_index_src(src);
            let index = csi::read(index_src)?;
            self.index = Some(Box::new(index));
        }

        let file = File::open(src)?;
        self.build_from_reader(file)
    }

    /// Builds a indexed SAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// use noodles_sam::indexed_reader::Builder;
    ///
    /// let index = csi::Index::default();
    /// let reader = Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<R>>
    where
        R: Read,
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
            build_index_src("sample.sam.gz"),
            PathBuf::from("sample.sam.gz.csi")
        );
    }
}
