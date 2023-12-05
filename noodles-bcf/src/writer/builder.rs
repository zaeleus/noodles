use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Writer;

/// A BCF writer builder.
pub struct Builder;

impl Builder {
    /// Builds a BCF writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::writer::Builder;
    /// let writer = Builder.build_from_path("out.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<bgzf::Writer<File>>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(Writer::new)
    }

    /// Builds a BCF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_bcf::writer::Builder;
    /// let writer = Builder.build_from_writer(io::sink())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> io::Result<Writer<bgzf::Writer<W>>>
    where
        W: Write,
    {
        Ok(Writer::new(writer))
    }
}
