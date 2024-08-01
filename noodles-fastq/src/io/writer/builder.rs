use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use super::Writer;

/// A FASTQ writer builder.
#[derive(Default)]
pub struct Builder;

impl Builder {
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
        Writer::new(Box::new(writer))
    }
}
