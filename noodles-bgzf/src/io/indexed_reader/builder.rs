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
    pub fn set_index(mut self, index: gzi::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed BGZF reader from a path.
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
