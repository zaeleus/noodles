use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use super::{DEFAULT_DEFINITION_SEPARATOR, Writer};

/// A FASTQ writer builder.
pub struct Builder {
    definition_separator: u8,
}

impl Builder {
    /// Sets the definition separator.
    ///
    /// By default, this is a space (` `).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::io::writer::Builder;
    /// let builder = Builder::default().set_definition_separator(b'\t');
    /// ```
    pub fn set_definition_separator(mut self, definition_separator: u8) -> Self {
        self.definition_separator = definition_separator;
        self
    }

    /// Builds a FASTQ writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fastq::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.fq")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<Box<dyn Write>>>
    where
        P: AsRef<Path>,
    {
        let writer = File::create(dst).map(BufWriter::new)?;
        Ok(self.build_from_writer(writer))
    }

    /// Builds a FASTQ writer from a writer.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_fastq::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(io::sink());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<Box<dyn Write>>
    where
        W: Write + 'static,
    {
        Writer {
            inner: Box::new(writer),
            definition_separator: self.definition_separator,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            definition_separator: DEFAULT_DEFINITION_SEPARATOR,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();
        assert_eq!(builder.definition_separator, b' ');
    }
}
