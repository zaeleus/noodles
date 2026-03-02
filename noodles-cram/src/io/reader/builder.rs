use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_fasta as fasta;

use super::Reader;
use crate::file_definition::Version;

/// A CRAM reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::reader::Builder;
    /// use noodles_fasta as fasta;
    ///
    /// let reference_sequence_repository = fasta::Repository::default();
    ///
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

    /// Builds a CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(|file| self.build_from_reader(file))
    }

    /// Builds a CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        Reader {
            inner: reader,
            reference_sequence_repository: self.reference_sequence_repository,
            version: Version::default(),
        }
    }
}
