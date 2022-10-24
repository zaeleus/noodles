use std::io::Write;

use super::Writer;

const DEFAULT_LINE_BASE_COUNT: usize = 80;

/// A FASTA writer builder.
pub struct Builder {
    line_base_count: usize,
}

impl Builder {
    /// Sets the number of bases per line.
    ///
    /// By default, this is set to 80.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let builder = fasta::writer::Builder::default().set_line_base_count(100);
    /// ```
    pub fn set_line_base_count(mut self, line_base_count: usize) -> Self {
        self.line_base_count = line_base_count;
        self
    }

    /// Builds a FASTA writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    /// let writer = fasta::writer::Builder::default().build_with_writer(io::sink());
    /// ```
    pub fn build_with_writer<W>(self, writer: W) -> Writer<W>
    where
        W: Write,
    {
        Writer {
            inner: writer,
            line_base_count: self.line_base_count,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            line_base_count: DEFAULT_LINE_BASE_COUNT,
        }
    }
}
