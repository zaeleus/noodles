use std::io::Write;

use noodles_fasta as fasta;

use super::{Options, Writer};
use crate::DataContainer;

/// A CRAM writer builder.
pub struct Builder<W> {
    inner: W,
    reference_sequences: Vec<fasta::Record>,
    options: Options,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(crate) fn new(inner: W, reference_sequences: Vec<fasta::Record>) -> Self {
        Self {
            inner,
            reference_sequences,
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

    /// Sets whether an external reference sequence is required to fully decode the record.
    ///
    /// If `false`, the sequence is embedded in the slice.
    ///
    /// The default is `true`.
    pub fn require_reference_sequence(mut self, value: bool) -> Self {
        self.options.require_reference_sequence = value;
        self
    }

    /// Builds a CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::builder(Vec::new(), Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        Writer {
            inner: self.inner,
            reference_sequences: self.reference_sequences,
            options: self.options,
            data_container_builder: DataContainer::builder(0),
            record_counter: 0,
        }
    }
}
