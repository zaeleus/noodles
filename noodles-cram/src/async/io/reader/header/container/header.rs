use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::r#async::io::reader::{
    num::{read_itf8, read_ltf8},
    CrcReader,
};

pub async fn read_header<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let mut crc_reader = CrcReader::new(reader);

    let length = crc_reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    // reference sequence ID
    read_itf8(&mut crc_reader).await?;

    // alignment start
    read_itf8(&mut crc_reader).await?;

    // alignment span
    read_itf8(&mut crc_reader).await?;

    // record count
    read_itf8(&mut crc_reader).await?;

    // record counter
    read_ltf8(&mut crc_reader).await?;

    // base count
    read_ltf8(&mut crc_reader).await?;

    // block count
    read_itf8(&mut crc_reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    read_landmarks(&mut crc_reader).await?;

    let actual_crc32 = crc_reader.crc().sum();

    let reader = crc_reader.into_inner();
    let expected_crc32 = reader.read_u32_le().await?;

    if actual_crc32 != expected_crc32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "container header checksum mismatch: expected {expected_crc32:08x}, got {actual_crc32:08x}"
            ),
        ));
    }

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
