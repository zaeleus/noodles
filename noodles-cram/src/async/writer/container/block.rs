use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::container::Block;

pub async fn write_block<W>(writer: &mut W, block: &Block) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::writer::container::write_block(&mut buf, block)?;
    writer.write_all(&buf).await?;
    Ok(())
}
