use std::io::Write;

use super::{CompressionLevel, MAX_BUF_SIZE, Writer};

/// A BGZF writer builder.
#[derive(Debug, Default)]
pub struct Builder {
    compression_level: CompressionLevel,
}

impl Builder {
    /// Sets a compression level.
    ///
    /// By default, the compression level is set to level 6.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, io::writer::CompressionLevel};
    ///
    /// let builder = bgzf::io::writer::Builder::default()
    ///     .set_compression_level(CompressionLevel::BEST);
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = compression_level;
        self
    }

    /// Builds a BGZF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::io::writer::Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<W>
    where
        W: Write,
    {
        Writer {
            inner: Some(writer),
            position: 0,
            staging_buf: Vec::with_capacity(MAX_BUF_SIZE),
            compression_buf: Vec::new(),
            compression_level: self.compression_level.into(),
        }
    }
}
