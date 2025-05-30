use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    r#async::io::reader::{
        CrcReader,
        num::{read_itf8, read_itf8_as, read_ltf8_as},
    },
    container::{Header, ReferenceSequenceContext},
};

pub async fn read_header<R>(reader: &mut R, header: &mut Header) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let mut crc_reader = CrcReader::new(reader);
    read_header_inner(&mut crc_reader, header).await
}

pub async fn read_header_inner<R>(
    reader: &mut CrcReader<R>,
    header: &mut Header,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::container::header::is_eof;

    let len = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence_id = read_itf8(reader).await?;
    let alignment_start = read_itf8(reader).await?;
    let alignment_span = read_itf8(reader).await?;

    header.reference_sequence_context = ReferenceSequenceContext::try_from((
        reference_sequence_id,
        alignment_start,
        alignment_span,
    ))?;

    header.record_count = read_itf8_as(reader).await?;
    header.record_counter = read_ltf8_as(reader).await?;
    header.base_count = read_ltf8_as(reader).await?;
    header.block_count = read_itf8_as(reader).await?;

    read_landmarks(reader, &mut header.landmarks).await?;

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

    if is_eof(
        len,
        reference_sequence_id,
        alignment_start,
        header.block_count,
        actual_crc32,
    ) {
        Ok(0)
    } else {
        Ok(len)
    }
}

async fn read_landmarks<R>(reader: &mut R, landmarks: &mut Vec<usize>) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    landmarks.clear();

    let len: usize = read_itf8_as(reader).await?;

    for _ in 0..len {
        let pos = read_itf8_as(reader).await?;
        landmarks.push(pos);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::container::ReferenceSequenceContext;

    #[tokio::test]
    async fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
        let src = [
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

        let mut actual = Header::default();
        let len = read_header(&mut &src[..], &mut actual).await?;

        let expected = Header {
            reference_sequence_context: ReferenceSequenceContext::some(
                2,
                Position::try_from(3)?,
                Position::try_from(7)?,
            ),
            record_count: 8,
            record_counter: 13,
            base_count: 21,
            block_count: 34,
            landmarks: vec![55, 89],
        };

        assert_eq!(len, 144);
        assert_eq!(actual, expected);

        Ok(())
    }
}
