use bytes::{Bytes, BytesMut};
use flate2::Compression;
use futures::SinkExt;
use tokio::io::AsyncWrite;
use tokio_util::codec::FramedWrite;

use super::{Deflater, Writer};
use crate::{block, r#async::BlockCodec, writer::BGZF_EOF};

/// An async BGZF writer builder.
#[derive(Debug)]
pub struct Builder<W> {
    inner: W,
    compression_level: Option<Compression>,
    worker_count: Option<usize>,
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
    /// use flate2::Compression;
    /// use noodles_bgzf as bgzf;
    ///
    /// let builder = bgzf::AsyncWriter::builder(Vec::new())
    ///     .set_compression_level(Compression::best());
    /// ```
    pub fn set_compression_level(&mut self, compression_level: Compression) -> &mut Self {
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
    /// use noodles_bgzf as bgzf;
    /// let builder = bgzf::AsyncWriter::builder(Vec::new()).set_worker_count(8);
    /// ```
    pub fn set_worker_count(&mut self, worker_count: usize) -> &mut Self {
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
        let worker_count = self.worker_count.unwrap_or_else(num_cpus::get);

        Writer {
            sink: Deflater::new(FramedWrite::new(self.inner, BlockCodec)).buffer(worker_count),
            buf: BytesMut::with_capacity(block::MAX_UNCOMPRESSED_DATA_LENGTH),
            eof_buf: Bytes::from_static(BGZF_EOF),
            compression_level,
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
