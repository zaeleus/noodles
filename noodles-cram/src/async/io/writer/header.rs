use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use noodles_sam as sam;

pub async fn write_header_container<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::header::write_header_container(&mut buf, header)?;
    writer.write_all(&buf).await?;
    Ok(())
}
