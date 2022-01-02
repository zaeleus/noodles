use noodles_bgzf::{self as bgzf, writer::CompressionLevel};
use tokio::io::AsyncWrite;

use super::Writer;

/// An async BAM writer builder.
pub struct Builder<W> {
    inner: W,
    compression_level: Option<CompressionLevel>,
    worker_count: Option<usize>,
}

impl<W> Builder<W>
where
    W: AsyncWrite + Unpin,
{
    pub(crate) fn new(inner: W) -> Self {
        Self {
            inner,
            compression_level: None,
            worker_count: None,
        }
    }

    /// Sets a compression level.
    ///
    /// By default, the compression level is set to level 6.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_bgzf::writer::CompressionLevel;
    ///
    /// let builder = bam::AsyncWriter::builder(Vec::new())
    ///     .set_compression_level(CompressionLevel::best());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = Some(compression_level);
        self
    }

    /// Sets a worker count.
    ///
    /// By default, the worker count is set to the number of available logical CPUs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let builder = bam::AsyncWriter::builder(Vec::new()).set_worker_count(8);
    /// ```
    pub fn set_worker_count(mut self, worker_count: usize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::AsyncWriter::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<bgzf::AsyncWriter<W>> {
        let mut builder = bgzf::AsyncWriter::builder(self.inner);

        if let Some(compression_level) = self.compression_level {
            builder = builder.set_compression_level(compression_level);
        }

        if let Some(worker_count) = self.worker_count {
            builder = builder.set_worker_count(worker_count);
        }

        Writer::from(builder.build())
    }
}
