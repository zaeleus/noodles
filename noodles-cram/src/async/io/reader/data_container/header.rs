use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    data_container::Header,
    r#async::io::reader::{
        num::{read_itf8, read_ltf8},
        CrcReader,
    },
};

pub async fn read_header<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::data_container::header::{build_reference_sequence_context, is_eof};

    let mut crc_reader = CrcReader::new(reader);

    let length = crc_reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence_id = read_itf8(&mut crc_reader).await?;
    let alignment_start = read_itf8(&mut crc_reader).await?;
    let alignment_span = read_itf8(&mut crc_reader).await?;

    let number_of_records = read_itf8(&mut crc_reader).await?;

    let record_counter = read_ltf8(&mut crc_reader).await.and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let bases = read_ltf8(&mut crc_reader).await.and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let number_of_blocks = read_itf8(&mut crc_reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let landmarks = read_landmarks(&mut crc_reader).await?;

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

    if is_eof(
        length,
        reference_sequence_id,
        alignment_start,
        number_of_blocks,
        actual_crc32,
    ) {
        return Ok(None);
    }

    let reference_sequence_context =
        build_reference_sequence_context(reference_sequence_id, alignment_start, alignment_span)?;

    let header = Header::builder()
        .set_length(length)
        .set_reference_sequence_context(reference_sequence_context)
        .set_record_count(number_of_records)
        .set_record_counter(record_counter)
        .set_base_count(bases)
        .set_block_count(number_of_blocks)
        .set_landmarks(landmarks)
        .build();

    Ok(Some(header))
}

async fn read_landmarks<R>(reader: &mut R) -> io::Result<Vec<usize>>
where
    R: AsyncRead + Unpin,
{
    let len = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut landmarks = Vec::with_capacity(len);

    for _ in 0..len {
        let pos = read_itf8(reader).await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        landmarks.push(pos);
    }

    Ok(landmarks)
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::data_container::ReferenceSequenceContext;

    #[tokio::test]
    async fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
            0x90, 0x00, 0x00, 0x00, // length = 144 bytes
            0x02, // reference sequence ID = 2
            0x03, // starting position on the reference = 3
            0x05, // alignment span = 5
            0x08, // number of records = 8
            0x0d, // record counter = 13
            0x15, // bases = 21
            0x22, // number of blocks = 34
            0x02, // landmark count = 2
            0x37, // landmarks[0] = 55
            0x59, // landmarks[1] = 89
            0x21, 0xf7, 0x9c, 0xed, // CRC32
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let expected = Header::builder()
            .set_length(144)
            .set_reference_sequence_context(ReferenceSequenceContext::some(
                2,
                Position::try_from(3)?,
                Position::try_from(7)?,
            ))
            .set_record_count(8)
            .set_record_counter(13)
            .set_base_count(21)
            .set_block_count(34)
            .set_landmarks(vec![55, 89])
            .build();

        assert_eq!(actual, Some(expected));

        Ok(())
    }
}
