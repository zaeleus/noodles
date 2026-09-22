use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_csi::BinningIndex;

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
    /// use noodles_sam::io::indexed_reader::Builder;
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
    /// use noodles_sam::io::indexed_reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.sam.gz")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let index = match self.index {
            Some(index) => index,
            None => crate::fs::read_associated_index(src).map(Box::new)?,
        };

        let file = File::open(src)?;

        Ok(IndexedReader::new(file, index))
    }

    /// Builds a indexed SAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// use noodles_sam::io::indexed_reader::Builder;
    ///
    /// let index = csi::Index::default();
    /// let reader = Builder::default().set_index(index).build_from_reader(io::empty())?;
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
