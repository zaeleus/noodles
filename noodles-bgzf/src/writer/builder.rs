use std::io::Write;

use super::{CompressionLevel, Writer, MAX_BUF_SIZE};

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
    /// use noodles_bgzf::{self as bgzf, writer::CompressionLevel};
    ///
    /// let builder = bgzf::writer::Builder::default()
    ///     .set_compression_level(CompressionLevel::best());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = compression_level;
        self
    }

    /// Builds a BGZF writer from a writer..
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::writer::Builder::default().build_with_writer(io::sink());
    /// ```
    pub fn build_with_writer<W>(self, writer: W) -> Writer<W>
    where
        W: Write,
    {
        Writer {
            inner: Some(writer),
            position: 0,
            buf: Vec::with_capacity(MAX_BUF_SIZE),
            compression_level: self.compression_level.into(),
        }
    }
}
