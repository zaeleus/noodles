use std::io::Write;

use super::Writer;

const DEFAULT_LINE_BASE_COUNT: usize = 80;

/// A FASTA writer builder.
pub struct Builder<W> {
    inner: W,
    line_base_count: usize,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(super) fn new(inner: W) -> Self {
        Builder {
            inner,
            line_base_count: DEFAULT_LINE_BASE_COUNT,
        }
    }

    /// Sets the number of bases per line.
    ///
    /// By default, this is set to 80.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let builder = fasta::Writer::builder(Vec::new()).set_line_base_count(100);
    /// ```
    pub fn set_line_base_count(mut self, line_base_count: usize) -> Self {
        self.line_base_count = line_base_count;
        self
    }

    /// Builds a FASTA writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let writer = fasta::Writer::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        Writer {
            inner: self.inner,
            line_base_count: self.line_base_count,
        }
    }
}
