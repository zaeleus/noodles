use std::{
    fs::File,
    io::{self, Write},
    num::NonZero,
    path::Path,
};

use super::Writer;

pub(crate) const DEFAULT_LINE_BASE_COUNT: NonZero<usize> = NonZero::new(80).unwrap();

/// A FASTA writer builder.
pub struct Builder {
    line_base_count: NonZero<usize>,
}

impl Builder {
    /// Sets the number of bases per line.
    ///
    /// By default, this is set to 80.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::num::NonZero;
    /// use noodles_fasta::io::writer::Builder;
    /// let builder = Builder::default().set_line_base_count(const { NonZero::new(100).unwrap() });
    /// ```
    pub fn set_line_base_count(mut self, line_base_count: NonZero<usize>) -> Self {
        self.line_base_count = line_base_count;
        self
    }

    /// Builds a FASTA writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fasta::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.fa")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(|file| self.build_from_writer(file))
    }

    /// Builds a FASTA writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<W>
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
