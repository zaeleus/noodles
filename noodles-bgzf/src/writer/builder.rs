use std::io::Write;

use super::{CompressionLevel, Writer, MAX_BUF_SIZE};

/// A BGZF writer builder.
#[derive(Debug)]
pub struct Builder<W> {
    inner: W,
    compression_level: Option<CompressionLevel>,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(crate) fn new(inner: W) -> Self {
        Self {
            inner,
            compression_level: None,
        }
    }

    /// Sets a compression level.
    ///
    /// By default, the compression level is set to level 6.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, writer::CompressionLevel};
    ///
    /// let builder = bgzf::Writer::builder(Vec::new())
    ///     .set_compression_level(CompressionLevel::best());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = Some(compression_level);
        self
    }

    /// Builds an async BGZF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        let compression_level = self.compression_level.unwrap_or_default();

        Writer {
            inner: Some(self.inner),
            position: 0,
            buf: Vec::with_capacity(MAX_BUF_SIZE),
            compression_level: compression_level.into(),
        }
    }
}
