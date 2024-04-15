use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub async fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::io::writer::container::EOF;
    writer.write_all(&EOF).await
}
