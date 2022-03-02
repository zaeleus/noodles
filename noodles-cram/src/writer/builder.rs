use std::io::Write;

use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{Options, Writer};
use crate::DataContainer;

/// A CRAM writer builder.
pub struct Builder<'a, W, A> {
    inner: W,
    reference_sequence_repository: fasta::Repository<A>,
    header: &'a sam::Header,
    options: Options,
}

impl<'a, W, A> Builder<'a, W, A>
where
    W: Write,
    A: fasta::repository::Adapter,
{
    pub(crate) fn new(
        inner: W,
        reference_sequence_repository: fasta::Repository<A>,
        header: &'a sam::Header,
    ) -> Self {
        Self {
            inner,
            reference_sequence_repository,
            header,
            options: Options::default(),
        }
    }

    /// Sets whether to preserve read names.
    ///
    /// If `false`, read names are discared.
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
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::new(Vec::new());
    /// let header = sam::Header::default();
    /// let writer = cram::Writer::builder(Vec::new(), repository, &header).build();
    /// ```
    pub fn build(self) -> Writer<'a, W, A> {
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
