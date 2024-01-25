use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Writer;
use crate::io::CompressionMethod;

/// A BCF writer builder.
#[derive(Debug, Default)]
pub struct Builder {
    compression_method: Option<CompressionMethod>,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::io::{writer::Builder, CompressionMethod};
    /// let builder = Builder::default().set_compression_method(CompressionMethod::Bgzf);
    /// ```
    pub fn set_compression_method(mut self, compression_method: CompressionMethod) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Builds a BCF writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<Box<dyn Write>>>
    where
        P: AsRef<Path>,
    {
        let file = File::create(dst)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds a BCF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<'w, W>(self, writer: W) -> Writer<Box<dyn Write + 'w>>
    where
        W: Write + 'w,
    {
        let inner: Box<dyn Write> = match self.compression_method {
            Some(CompressionMethod::Bgzf) | None => Box::new(bgzf::Writer::new(writer)),
            Some(CompressionMethod::None) => Box::new(BufWriter::new(writer)),
        };

        Writer::from(inner)
    }
}
