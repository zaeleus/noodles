use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_fasta as fasta;

use super::IndexedReader;
use crate::crai;

/// An indexed CRAM reader builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    index: Option<crai::Index>,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::indexed_reader::Builder;
    /// use noodles_fasta as fasta;
    ///
    /// let reference_sequence_repository = fasta::Repository::default();
    /// let builder = Builder::default()
    ///     .set_reference_sequence_repository(reference_sequence_repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets an index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{crai, io::indexed_reader::Builder};
    /// let index = crai::Index::default();
    /// let builder = Builder::default().set_index(index);
    /// ```
    pub fn set_index(mut self, index: crai::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed CRAM reader from a path.
    ///
    /// If no index is set, this will attempt to read an associated index at `<src>.crai`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_cram::io::indexed_reader::Builder;
    /// let builder = Builder::default().build_from_path("sample.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.index.is_none() {
            self.index = crate::fs::read_associated_index(src).map(Some)?;
        }

        let file = File::open(src)?;
        self.build_from_reader(file)
    }

    /// Builds an indexed CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::{crai, io::indexed_reader::Builder};
    ///
    /// let index = crai::Index::default();
    /// let builder = Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(io::empty())?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<R>>
    where
        R: Read,
    {
        let inner = crate::io::reader::Builder::default()
            .set_reference_sequence_repository(self.reference_sequence_repository)
            .build_from_reader(reader);

        let index = self
            .index
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index"))?;

        Ok(IndexedReader { inner, index })
    }
}
