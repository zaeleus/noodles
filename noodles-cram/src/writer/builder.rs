use std::io::Write;

use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{Options, Writer};
use crate::DataContainer;

/// A CRAM writer builder.
pub struct Builder<'a, W> {
    inner: W,
    reference_sequence_repository: fasta::Repository,
    header: &'a sam::Header,
    options: Options,
}

impl<'a, W> Builder<'a, W>
where
    W: Write,
{
    pub(crate) fn new(inner: W, header: &'a sam::Header) -> Self {
        Self {
            inner,
            reference_sequence_repository: fasta::Repository::default(),
            header,
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
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::default();
    /// let writer = cram::Writer::builder(Vec::new(), &header).build();
    /// ```
    pub fn build(self) -> Writer<'a, W> {
        Writer {
            inner: self.inner,
            reference_sequence_repository: self.reference_sequence_repository,
            header: self.header,
            options: self.options,
            data_container_builder: DataContainer::builder(0),
            record_counter: 0,
        }
    }
}
