use std::num::NonZeroUsize;

use bytes::{Bytes, BytesMut};
use futures::SinkExt;
use tokio::io::AsyncWrite;
use tokio_util::codec::FramedWrite;

use super::{Deflater, Writer};
use crate::{
    r#async::BlockCodec,
    writer::{CompressionLevel, BGZF_EOF, MAX_BUF_SIZE},
};

/// An async BGZF writer builder.
#[derive(Debug)]
pub struct Builder<W> {
    inner: W,
    compression_level: Option<CompressionLevel>,
    worker_count: Option<NonZeroUsize>,
}

impl<W> Builder<W>
where
    W: AsyncWrite,
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
    /// use noodles_bgzf::{self as bgzf, writer::CompressionLevel};
    ///
    /// let builder = bgzf::AsyncWriter::builder(Vec::new())
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_bgzf as bgzf;
    ///
    /// let worker_count = NonZeroUsize::try_from(1)?;
    /// let builder = bgzf::AsyncWriter::builder(Vec::new()).set_worker_count(worker_count);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BGZF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::AsyncWriter::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        let compression_level = self.compression_level.unwrap_or_default();

        let worker_count = self.worker_count.unwrap_or_else(|| {
            // SAFETY: `num_cpus::get` is guaranteed to be non-zero.
            NonZeroUsize::new(num_cpus::get()).unwrap()
        });

        Writer {
            sink: Deflater::new(FramedWrite::new(self.inner, BlockCodec))
                .buffer(worker_count.get()),
            buf: BytesMut::with_capacity(MAX_BUF_SIZE),
            eof_buf: Bytes::from_static(BGZF_EOF),
            compression_level: compression_level.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let builder = Builder::new(Vec::new());
        assert!(builder.compression_level.is_none());
        assert!(builder.worker_count.is_none());
    }
}
