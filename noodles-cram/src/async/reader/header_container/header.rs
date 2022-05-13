use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::r#async::reader::num::{read_itf8, read_ltf8};

pub async fn read_header<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let length = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    // reference sequence ID
    read_itf8(reader).await?;

    // alignment start
    read_itf8(reader).await?;

    // alignment span
    read_itf8(reader).await?;

    // record count
    read_itf8(reader).await?;

    // record counter
    read_ltf8(reader).await?;

    // base count
    read_ltf8(reader).await?;

    // block count
    read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    read_landmarks(reader).await?;

    Ok(length)
}

async fn read_landmarks<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let len = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    for _ in 0..len {
        read_itf8(reader).await?;
    }

    Ok(())
}
