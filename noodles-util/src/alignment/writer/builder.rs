use std::io::Write;

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::Format;

/// An alignment writer builder.
pub struct Builder<W> {
    inner: W,
    format: Format,
    reference_sequence_repository: fasta::Repository,
}

impl<W> Builder<W>
where
    W: Write + 'static,
{
    pub(super) fn new(inner: W) -> Self {
        Self {
            inner,
            format: Format::Sam,
            reference_sequence_repository: fasta::Repository::default(),
        }
    }

    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::Writer::builder(io::sink()).set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::Writer::builder(io::sink())
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an alignment writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let writer = alignment::Writer::builder(io::sink())
    ///     .set_format(Format::Sam)
    ///     .build();
    /// ```
    pub fn build(self) -> Writer {
        let inner: Box<dyn sam::AlignmentWriter> = match self.format {
            Format::Sam => Box::new(sam::Writer::new(self.inner)),
            Format::Bam => Box::new(bam::Writer::new(self.inner)),
            Format::Cram => Box::new(
                cram::Writer::builder(self.inner)
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .build(),
            ),
        };

        Writer { inner }
    }
}
