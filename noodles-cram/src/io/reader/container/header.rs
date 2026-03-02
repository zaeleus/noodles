use std::io::{self, Read};

use flate2::CrcReader;

use crate::{
    container::{Header, ReferenceSequenceContext},
    file_definition::Version,
    io::reader::num::{
        read_header_int, read_i32_le, read_long_as, read_position, read_u32_le, read_uint7_as,
        read_unsigned_int_as,
    },
};

// § 9 "End of file container" (2022-04-12)
const EOF_LENGTH: usize = 15;
// CRAM 2.x EOF has length 11 (no CRC32 in the block: 15 - 4 = 11)
const EOF_LENGTH_V2: usize = 11;
const EOF_REFERENCE_SEQUENCE_ID: i32 = -1;
// CRAM 4.0 EOF: htslib encodes ref_seq_id = -1 as zigzag(1) in the hardcoded EOF bytes,
// but `read_header_int` decodes it as unsigned (1), not as signed (-1).
const EOF_REFERENCE_SEQUENCE_ID_V4: i32 = 1;
const EOF_ALIGNMENT_START: i64 = 4_542_278;
const EOF_BLOCK_COUNT: usize = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;
// CRAM 4.0 EOF has a different CRC32 because the header uses VLQ encoding
const EOF_CRC32_V4: u32 = 0xae_f7_8f_52;

pub(super) fn read_header<R>(
    reader: &mut R,
    header: &mut Header,
    version: Version,
) -> io::Result<usize>
where
    R: Read,
{
    let mut crc_reader = CrcReader::new(reader);

    match read_header_inner(&mut crc_reader, header, version) {
        Ok(len) => Ok(len),
        // An `UnexpectedEof` during container header reading means there is no more data.
        // CRAM 2.0 has no EOF marker (introduced in 2.1), but this can also occur when
        // seeking past the EOF container in other versions (e.g., `query_unmapped()`
        // seeking to `End(0)`).
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

fn read_header_inner<R>(
    reader: &mut CrcReader<R>,
    header: &mut Header,
    version: Version,
) -> io::Result<usize>
where
    R: Read,
{
    let len = if version.uses_vlq() {
        read_uint7_as(reader)?
    } else {
        read_i32_le(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?
    };

    let reference_sequence_id = read_header_int(reader, version)?;

    let alignment_start = read_position(reader, version)?;
    let alignment_span = read_position(reader, version)?;

    header.record_count = read_unsigned_int_as(reader, version)?;
    header.record_counter = read_long_as(reader, version)?;
    header.base_count = read_long_as(reader, version)?;
    header.block_count = read_unsigned_int_as(reader, version)?;

    read_landmarks(reader, &mut header.landmarks, version)?;

    if version.has_crc32() {
        let actual_crc32 = reader.crc().sum();
        let expected_crc32 = read_u32_le(reader.get_mut())?;

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
            expected_crc32,
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

fn read_landmarks<R>(reader: &mut R, landmarks: &mut Vec<usize>, version: Version) -> io::Result<()>
where
    R: Read,
{
    landmarks.clear();

    let n: usize = read_unsigned_int_as(reader, version)?;

    for _ in 0..n {
        let pos = read_unsigned_int_as(reader, version)?;
        landmarks.push(pos);
    }

    Ok(())
}

pub(crate) fn is_eof(
    length: usize,
    reference_sequence_id: i32,
    alignment_start: i64,
    block_count: usize,
    crc32: u32,
) -> bool {
    let is_v3_eof = length == EOF_LENGTH
        && reference_sequence_id == EOF_REFERENCE_SEQUENCE_ID
        && alignment_start == EOF_ALIGNMENT_START
        && block_count == EOF_BLOCK_COUNT
        && crc32 == EOF_CRC32;

    let is_v4_eof = length == EOF_LENGTH
        && reference_sequence_id == EOF_REFERENCE_SEQUENCE_ID_V4
        && alignment_start == EOF_ALIGNMENT_START
        && block_count == EOF_BLOCK_COUNT
        && crc32 == EOF_CRC32_V4;

    is_v3_eof || is_v4_eof
}

pub(crate) fn is_eof_v2(
    length: usize,
    reference_sequence_id: i32,
    alignment_start: i64,
    block_count: usize,
) -> bool {
    length == EOF_LENGTH_V2
        && reference_sequence_id == EOF_REFERENCE_SEQUENCE_ID
        && alignment_start == EOF_ALIGNMENT_START
        && block_count == EOF_BLOCK_COUNT
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::file_definition::Version;

    #[test]
    fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
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
        let len = read_header(&mut &src[..], &mut actual, Version::V3_0)?;

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

    #[test]
    fn test_read_header_with_eof() -> io::Result<()> {
        let src = [
            0x0f, 0x00, 0x00, 0x00, // length = 15 bytes
            0xff, 0xff, 0xff, 0xff, 0x0f, // reference sequence ID = None (-1)
            0xe0, 0x45, 0x4f, 0x46, // starting position on the reference = 4542278
            0x00, // alignment span = 0
            0x00, // number of records = 0
            0x00, // record counter = 0
            0x00, // bases = 0
            0x01, // number of blocks = 1
            0x00, // landmark count = 0
            0x05, 0xbd, 0xd9, 0x4f, // CRC32
        ];

        let mut header = Header::default();
        let len = read_header(&mut &src[..], &mut header, Version::V3_0)?;

        assert_eq!(len, 0);

        Ok(())
    }

    #[test]
    fn test_read_header_with_a_checksum_mismatch() {
        // EOF container header
        let src = [
            0x0f, 0x00, 0x00, 0x00, // length = 15 bytes
            0xff, 0xff, 0xff, 0xff, 0x0f, // reference sequence ID = None (-1)
            0xe0, 0x45, 0x4f, 0x46, // starting position on the reference = 4542278
            0x00, // alignment span = 0
            0x00, // number of records = 0
            0x00, // record counter = 0
            0x00, // bases = 0
            0x01, // number of blocks = 1
            0x00, // landmark count = 0
            0x00, 0x00, 0x00, 0x00, // CRC32 (invalid)
        ];

        let mut header = Header::default();

        assert!(matches!(
            read_header(&mut &src[..], &mut header, Version::V3_0),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }
}
