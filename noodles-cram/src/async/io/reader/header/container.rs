mod block;
mod header;

use tokio::io::{self, AsyncRead, AsyncReadExt, Take};

use self::block::read_block;
pub(super) use self::header::read_header;

pub(super) struct Reader<R> {
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

    pub(super) async fn raw_sam_header_reader(
        &mut self,
    ) -> io::Result<impl AsyncRead + Unpin + '_> {
        let mut reader = read_block(&mut self.inner).await?;

        let len = reader.read_i32_le().await.and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        Ok(reader.take(len))
    }
}
