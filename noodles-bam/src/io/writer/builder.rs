use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Writer;

/// A BAM writer builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a BAM writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bam as bam;
    /// let writer = bam::io::writer::Builder::default().build_from_path("out.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<bgzf::io::Writer<File>>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(Writer::new)
    }

    /// Builds a BAM writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let writer = bam::io::writer::Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<bgzf::io::Writer<W>>
    where
        W: Write,
    {
        Writer::new(writer)
    }
}
