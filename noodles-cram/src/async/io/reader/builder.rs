use std::path::Path;

use noodles_fasta as fasta;
use tokio::{
    fs::File,
    io::{self, AsyncRead},
};

use super::Reader;
use crate::file_definition::Version;

/// An async CRAM reader builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::reader::Builder;
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

    /// Builds an async CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram::r#async::io::reader::Builder;
    /// let _reader = Builder::default().build_from_path("sample.cram").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src)
            .await
            .map(|file| self.build_from_reader(file))
    }

    /// Builds an async CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::reader::Builder;
    /// use tokio::io;
    /// let _reader = Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: AsyncRead + Unpin,
    {
        Reader {
            inner: reader,
            reference_sequence_repository: self.reference_sequence_repository,
            version: Version::default(),
        }
    }
}
