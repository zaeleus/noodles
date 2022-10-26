mod block;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub use self::block::write_block;

pub async fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::writer::container::EOF;
    writer.write_all(&EOF).await
}
