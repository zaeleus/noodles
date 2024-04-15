use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::DataContainer;

pub async fn write_data_container<W>(
    writer: &mut W,
    data_container: &DataContainer,
    base_count: u64,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::data_container::write_data_container(&mut buf, data_container, base_count)?;
    writer.write_all(&buf).await?;
    Ok(())
}
