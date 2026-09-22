use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use super::IndexedReader;
use crate::{gzi, io::reader};

/// An indexed BGZF reader builder.
#[derive(Default)]
pub struct Builder {
    reader_builder: reader::Builder,
    index: Option<gzi::Index>,
}

impl Builder {
    /// Sets a GZ index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, gzi};
    /// let index = gzi::Index::default();
    /// let builder = bgzf::io::indexed_reader::Builder::default().set_index(index);
    /// ```
    pub fn set_index(mut self, index: gzi::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed BGZF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::indexed_reader::Builder::default().build_from_path("src.gz")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let index = match self.index {
            Some(index) => index,
            None => crate::fs::read_associated_index(src)?,
        };

        let inner = self.reader_builder.build_from_path(src)?;

        Ok(IndexedReader { inner, index })
    }

    /// Builds a indexed BGZF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_bgzf::{self as bgzf, gzi};
    ///
    /// let index = gzi::Index::default();
    /// let reader = bgzf::io::indexed_reader::Builder::default()
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

        let inner = self.reader_builder.build_from_reader(reader);

        Ok(IndexedReader { inner, index })
    }
}
