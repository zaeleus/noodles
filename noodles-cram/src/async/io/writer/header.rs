use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use noodles_sam as sam;

use crate::FileDefinition;

pub(super) async fn write_file_definition<W>(
    writer: &mut W,
    file_definition: &FileDefinition,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::header::write_file_definition(&mut buf, file_definition)?;
    writer.write_all(&buf).await?;
    Ok(())
}

pub(super) async fn write_file_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::header::write_file_header(&mut buf, header)?;
    writer.write_all(&buf).await?;
    Ok(())
}
