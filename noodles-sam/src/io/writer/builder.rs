use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Writer;
use crate::io::CompressionMethod;

/// A SAM writer builder.
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
    /// use noodles_sam::io::{writer::Builder, CompressionMethod};
    /// let builder = Builder::default().set_compression_method(CompressionMethod::Bgzf);
    /// ```
    pub fn set_compression_method(mut self, compression_method: CompressionMethod) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Builds a SAM writer from a path.
    ///
    /// If the compression method is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_sam::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("sample.sam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, dst: P) -> io::Result<Writer<Box<dyn Write>>>
    where
        P: AsRef<Path>,
    {
        let dst = dst.as_ref();

        if self.compression_method.is_none() {
            self.compression_method = match dst.extension().and_then(|ext| ext.to_str()) {
                Some("gz" | "bgz") => Some(CompressionMethod::Bgzf),
                _ => Some(CompressionMethod::None),
            };
        }

        let file = File::create(dst)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds a SAM writer from a writer.
    ///
    /// If the compression method is not set, no compression is used.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<'w, W>(self, writer: W) -> Writer<Box<dyn Write + 'w>>
    where
        W: Write + 'w,
    {
        let inner: Box<dyn Write> = match self.compression_method {
            Some(CompressionMethod::Bgzf) => Box::new(bgzf::io::Writer::new(writer)),
            Some(CompressionMethod::None) | None => Box::new(BufWriter::new(writer)),
        };

        Writer::new(inner)
    }
}
