use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::Container;

pub async fn write_container<W>(
    writer: &mut W,
    container: &Container,
    base_count: u64,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::container::write_container(&mut buf, container, base_count)?;
    writer.write_all(&buf).await?;
    Ok(())
}

pub async fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::io::writer::container::EOF;
    writer.write_all(&EOF).await
}
