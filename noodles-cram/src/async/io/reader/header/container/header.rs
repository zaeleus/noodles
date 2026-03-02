use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    r#async::io::reader::{
        CrcReader,
        num::{read_header_int, read_long_as, read_position, read_uint7, read_unsigned_int_as},
    },
    file_definition::Version,
};

pub(crate) async fn read_header<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    if version.has_crc32() {
        let mut crc_reader = CrcReader::new(reader);
        read_header_with_crc32(&mut crc_reader, version).await
    } else {
        read_header_without_crc32(reader, version).await
    }
}

async fn read_length<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    if version.uses_vlq() {
        read_uint7(reader).await.map(u64::from)
    } else {
        reader.read_i32_le().await.and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}

async fn read_header_with_crc32<R>(reader: &mut CrcReader<R>, version: Version) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    let length = read_length(reader, version).await?;

    let _reference_sequence_id = read_header_int(reader, version).await?;
    let _alignment_start = read_position(reader, version).await?;
    let _alignment_span = read_position(reader, version).await?;
    let _record_count: usize = read_unsigned_int_as(reader, version).await?;
    let _record_counter: usize = read_long_as(reader, version).await?;
    let _base_count: usize = read_long_as(reader, version).await?;
    let _block_count: usize = read_unsigned_int_as(reader, version).await?;
    read_landmarks(reader, version).await?;

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

async fn read_header_without_crc32<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    let length = read_length(reader, version).await?;

    let _reference_sequence_id = read_header_int(reader, version).await?;
    let _alignment_start = read_position(reader, version).await?;
    let _alignment_span = read_position(reader, version).await?;
    let _record_count: usize = read_unsigned_int_as(reader, version).await?;
    let _record_counter: usize = read_long_as(reader, version).await?;
    let _base_count: usize = read_long_as(reader, version).await?;
    let _block_count: usize = read_unsigned_int_as(reader, version).await?;
    read_landmarks(reader, version).await?;

    Ok(length)
}

async fn read_landmarks<R>(reader: &mut R, version: Version) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let len: usize = read_unsigned_int_as(reader, version).await?;

    for _ in 0..len {
        let _landmark: usize = read_unsigned_int_as(reader, version).await?;
    }

    Ok(())
}
