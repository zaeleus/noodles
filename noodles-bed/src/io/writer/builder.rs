use std::{
    fs::File,
    io::BufWriter,
    io::{self, Write},
    path::Path,
};

use super::Writer;

/// A BED writer builder.
#[derive(Default)]
pub struct Builder<const N: usize>;

impl<const N: usize> Builder<N> {
    /// Builds a BED writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bed::io::writer::Builder;
    /// let writer = Builder::<3>.build_from_path("out.bed")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<N, BufWriter<File>>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(BufWriter::new).map(Writer::new)
    }

    /// Builds a BED writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed::io::writer::Builder;
    /// let writer = Builder::<3>.build_from_writer(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<N, BufWriter<W>>
    where
        W: Write,
    {
        Writer::new(BufWriter::new(writer))
    }
}
