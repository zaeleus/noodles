//! Async CRAM header container reader.

mod block;
mod header;
pub mod sam_header;

use tokio::io::{self, AsyncRead, AsyncReadExt, Take};

pub use self::block::Decoder;
use self::block::read_block;
pub(super) use self::header::read_header;
use crate::file_definition::Version;

/// An async CRAM header container reader.
pub struct Reader<R> {
    inner: Take<R>,
    version: Version,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    pub(super) fn new(inner: R, len: u64, version: Version) -> Self {
        Self {
            inner: inner.take(len),
            version,
        }
    }

    /// Returns a raw SAM header reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header text, e.g., using
    /// [`sam_header::Reader::discard_to_end`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::{fs::File, io::AsyncReadExt};
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number().await?;
    /// header_reader.read_format_version().await?;
    /// header_reader.read_file_id().await?;
    ///
    /// let mut container_reader = header_reader.container_reader().await?;
    ///
    /// let buf = {
    ///     let mut buf = Vec::new();
    ///
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader().await?;
    ///     raw_sam_header_reader.read_to_end(&mut buf).await?;
    ///     raw_sam_header_reader.discard_to_end().await?;
    ///
    ///     buf
    /// };
    ///
    /// container_reader.discard_to_end().await?;
    /// # Ok::<_, std::io::Error>(())
    /// # }
    /// ```
    pub async fn raw_sam_header_reader(
        &mut self,
    ) -> io::Result<sam_header::Reader<Decoder<&mut Take<R>>>> {
        let mut reader = read_block(&mut self.inner, self.version).await?;

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
