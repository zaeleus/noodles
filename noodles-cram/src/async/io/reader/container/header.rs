use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    r#async::io::reader::{
        CrcReader,
        num::{read_header_int, read_long_as, read_position, read_uint7_as, read_unsigned_int_as},
    },
    container::{Header, ReferenceSequenceContext},
    file_definition::Version,
};

pub(super) async fn read_header<R>(
    reader: &mut R,
    header: &mut Header,
    version: Version,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let mut crc_reader = CrcReader::new(reader);

    match read_header_inner(&mut crc_reader, header, version).await {
        Ok(len) => Ok(len),
        // An `UnexpectedEof` during container header reading means there is no more data.
        // CRAM 2.0 has no EOF marker (introduced in 2.1), but this can also occur when
        // seeking past the EOF container in other versions (e.g., `query_unmapped()`
        // seeking to `End(0)`).
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

async fn read_header_inner<R>(
    reader: &mut CrcReader<R>,
    header: &mut Header,
    version: Version,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::container::header::{is_eof, is_eof_v2};

    let len = if version.uses_vlq() {
        read_uint7_as(reader).await?
    } else {
        reader.read_i32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?
    };

    let reference_sequence_id = read_header_int(reader, version).await?;
    let alignment_start = read_position(reader, version).await?;
    let alignment_span = read_position(reader, version).await?;

    header.record_count = read_unsigned_int_as(reader, version).await?;
    header.record_counter = read_long_as(reader, version).await?;
    header.base_count = read_long_as(reader, version).await?;
    header.block_count = read_unsigned_int_as(reader, version).await?;

    read_landmarks(reader, &mut header.landmarks, version).await?;

    if version.has_crc32() {
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
            return Ok(0);
        }
    } else {
        // CRAM 2.x: check EOF without CRC32
        if is_eof_v2(
            len,
            reference_sequence_id,
            alignment_start,
            header.block_count,
        ) {
            return Ok(0);
        }
    }

    // Construct ReferenceSequenceContext after EOF check, since the EOF container
    // has field values that don't form a valid context.
    header.reference_sequence_context = ReferenceSequenceContext::try_from((
        reference_sequence_id,
        alignment_start,
        alignment_span,
    ))?;

    Ok(len)
}

async fn read_landmarks<R>(
    reader: &mut R,
    landmarks: &mut Vec<usize>,
    version: Version,
) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    landmarks.clear();

    let len: usize = read_unsigned_int_as(reader, version).await?;

    for _ in 0..len {
        let pos = read_unsigned_int_as(reader, version).await?;
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
        let len = read_header(&mut &src[..], &mut actual, Version::V3_0).await?;

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
