use std::io::Write;

use noodles_fasta as fasta;

use super::{Options, Writer};
use crate::DataContainer;

/// A CRAM writer builder.
pub struct Builder<W> {
    inner: W,
    reference_sequence_repository: fasta::Repository,
    options: Options,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(crate) fn new(inner: W) -> Self {
        Self {
            inner,
            reference_sequence_repository: fasta::Repository::default(),
            options: Options::default(),
        }
    }

    /// Sets the reference sequence repository.
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets whether to preserve read names.
    ///
    /// If `false`, read names are discarded.
    ///
    /// The default is `true`.
    pub fn preserve_read_names(mut self, value: bool) -> Self {
        self.options.preserve_read_names = value;
        self
    }

    /// Sets whether to encode alignment start positions as deltas.
    ///
    /// If `false`, record alignment start positions are written with their actual values.
    ///
    /// The default is `true`.
    pub fn encode_alignment_start_positions_as_deltas(mut self, value: bool) -> Self {
        self.options.encode_alignment_start_positions_as_deltas = value;
        self
    }

    /// Builds a CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        Writer {
            inner: self.inner,
            reference_sequence_repository: self.reference_sequence_repository,
            options: self.options,
            data_container_builder: DataContainer::builder(0),
            record_counter: 0,
        }
    }
}
