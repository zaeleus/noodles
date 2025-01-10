use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::r#async::io::reader::{
    num::{read_itf8, read_itf8_as, read_ltf8},
    CrcReader,
};

pub(crate) async fn read_header<R>(reader: &mut R) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    let mut crc_reader = CrcReader::new(reader);
    read_header_inner(&mut crc_reader).await
}

async fn read_header_inner<R>(reader: &mut CrcReader<R>) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    let length = reader.read_i32_le().await.and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let _reference_sequence_id = read_itf8(reader).await?;
    let _alignment_start = read_itf8(reader).await?;
    let _alignment_span = read_itf8(reader).await?;
    let _record_count = read_itf8(reader).await?;
    let _record_counter = read_ltf8(reader).await?;
    let _base_count = read_ltf8(reader).await?;
    let _block_count = read_itf8(reader).await?;
    read_landmarks(reader).await?;

    let actual_crc32 = reader.crc().sum();
    let expected_crc32 = reader.get_mut().read_u32_le().await?;

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
    let len: usize = read_itf8_as(reader).await?;

    for _ in 0..len {
        read_itf8(reader).await?;
    }

    Ok(())
}
