use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub async fn write_header<W>(writer: &mut W, len: usize) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::writer::header_container::write_header(&mut buf, len)?;
    writer.write_all(&buf).await?;
    Ok(())
}
