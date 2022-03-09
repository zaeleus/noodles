use noodles_core::Position;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    container::{self, ReferenceSequenceId},
    r#async::reader::num::{read_itf8, read_ltf8},
};

pub async fn read_header<R>(reader: &mut R) -> io::Result<container::Header>
where
    R: AsyncRead + Unpin,
{
    let length = reader.read_i32_le().await?;

    let reference_sequence_id = read_itf8(reader).await.and_then(|n| {
        ReferenceSequenceId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let starting_position_on_the_reference = read_itf8(reader)
        .await
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(Position::new)?;

    let alignment_span = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let number_of_records = read_itf8(reader).await?;
    let record_counter = read_ltf8(reader).await?;
    let bases = read_ltf8(reader).await?;

    let number_of_blocks = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let landmarks = read_landmarks(reader).await?;
    let crc32 = reader.read_u32_le().await?;

    let mut builder = container::Header::builder()
        .set_length(length)
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_span(alignment_span)
        .set_record_count(number_of_records)
        .set_record_counter(record_counter)
        .set_base_count(bases)
        .set_block_count(number_of_blocks)
        .set_landmarks(landmarks)
        .set_crc32(crc32);

    if let Some(position) = starting_position_on_the_reference {
        builder = builder.set_start_position(position);
    }

    Ok(builder.build())
}

async fn read_landmarks<R>(reader: &mut R) -> io::Result<Vec<i32>>
where
    R: AsyncRead + Unpin,
{
    let len = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;
    let mut landmarks = Vec::with_capacity(len);

    for _ in 0..len {
        let pos = read_itf8(reader).await?;
        landmarks.push(pos);
    }

    Ok(landmarks)
}

#[cfg(test)]
mod tests {
    use super::*;

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
            0xb4, 0x9f, 0x9c, 0xda, // CRC32
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let expected = container::Header::builder()
            .set_length(144)
            .set_reference_sequence_id(ReferenceSequenceId::try_from(2)?)
            .set_start_position(Position::try_from(3)?)
            .set_alignment_span(5)
            .set_record_count(8)
            .set_record_counter(13)
            .set_base_count(21)
            .set_block_count(34)
            .set_landmarks(vec![55, 89])
            .set_crc32(0xda9c9fb4)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
