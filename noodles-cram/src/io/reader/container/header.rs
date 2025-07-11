use std::io::{self, Read};

use flate2::CrcReader;

use crate::{
    container::{Header, ReferenceSequenceContext},
    io::reader::num::{read_i32_le, read_itf8, read_itf8_as, read_ltf8_as, read_u32_le},
};

// § 9 "End of file container" (2022-04-12)
const EOF_LENGTH: usize = 15;
const EOF_REFERENCE_SEQUENCE_ID: i32 = -1;
const EOF_ALIGNMENT_START: i32 = 4_542_278;
const EOF_BLOCK_COUNT: usize = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

pub fn read_header<R>(reader: &mut R, header: &mut Header) -> io::Result<usize>
where
    R: Read,
{
    let mut crc_reader = CrcReader::new(reader);
    read_header_inner(&mut crc_reader, header)
}

pub fn read_header_inner<R>(reader: &mut CrcReader<R>, header: &mut Header) -> io::Result<usize>
where
    R: Read,
{
    let len = read_i32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence_id = read_itf8(reader)?;
    let alignment_start = read_itf8(reader)?;
    let alignment_span = read_itf8(reader)?;

    header.reference_sequence_context = ReferenceSequenceContext::try_from((
        reference_sequence_id,
        alignment_start,
        alignment_span,
    ))?;

    header.record_count = read_itf8_as(reader)?;
    header.record_counter = read_ltf8_as(reader)?;
    header.base_count = read_ltf8_as(reader)?;
    header.block_count = read_itf8_as(reader)?;

    read_landmarks(reader, &mut header.landmarks)?;

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
        Ok(0)
    } else {
        Ok(len)
    }
}

fn read_landmarks<R>(reader: &mut R, landmarks: &mut Vec<usize>) -> io::Result<()>
where
    R: Read,
{
    landmarks.clear();

    let n: usize = read_itf8_as(reader)?;

    for _ in 0..n {
        let pos = read_itf8_as(reader)?;
        landmarks.push(pos);
    }

    Ok(())
}

pub(crate) fn is_eof(
    length: usize,
    reference_sequence_id: i32,
    alignment_start: i32,
    block_count: usize,
    crc32: u32,
) -> bool {
    length == EOF_LENGTH
        && reference_sequence_id == EOF_REFERENCE_SEQUENCE_ID
        && alignment_start == EOF_ALIGNMENT_START
        && block_count == EOF_BLOCK_COUNT
        && crc32 == EOF_CRC32
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

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
        let len = read_header(&mut &src[..], &mut actual)?;

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
        let len = read_header(&mut &src[..], &mut header)?;

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
            read_header(&mut &src[..], &mut header),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }
}
