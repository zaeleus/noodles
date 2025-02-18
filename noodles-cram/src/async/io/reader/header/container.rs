//! Async CRAM header container reader.

mod block;
mod header;
mod sam_header;

use tokio::io::{self, AsyncRead, AsyncReadExt, Take};

use self::block::read_block;
pub use self::block::Decoder;
pub(super) use self::header::read_header;

/// An async CRAM header container reader.
pub struct Reader<R> {
    inner: Take<R>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    pub(super) fn new(inner: R, len: u64) -> Self {
        Self {
            inner: inner.take(len),
        }
    }

    /// Returns a raw SAM header reader.
    pub async fn raw_sam_header_reader(
        &mut self,
    ) -> io::Result<sam_header::Reader<Decoder<&mut Take<R>>>> {
        let mut reader = read_block(&mut self.inner).await?;

        let len = reader.read_i32_le().await.and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        Ok(sam_header::Reader::new(reader, len))
    }

    /// Discards all input until EOF.
    pub async fn discard_to_end(&mut self) -> io::Result<u64> {
        io::copy(&mut self.inner, &mut io::sink()).await
    }
}
