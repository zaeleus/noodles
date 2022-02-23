use std::io::Write;

use noodles_fasta as fasta;

use super::Writer;
use crate::DataContainer;

/// A CRAM writer builder.
pub struct Builder<W> {
    inner: W,
    reference_sequences: Vec<fasta::Record>,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(crate) fn new(inner: W, reference_sequences: Vec<fasta::Record>) -> Self {
        Self {
            inner,
            reference_sequences,
        }
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
            data_container_builder: DataContainer::builder(0),
            record_counter: 0,
        }
    }
}
