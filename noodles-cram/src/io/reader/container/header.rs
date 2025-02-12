use std::{
    io::{self, Read},
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use flate2::CrcReader;
use noodles_core::Position;

use crate::{
    container::{Header, ReferenceSequenceContext},
    io::reader::num::{read_itf8, read_itf8_as, read_ltf8},
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
    let len = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence_id = read_itf8(reader)?;
    let alignment_start = read_itf8(reader)?;
    let alignment_span = read_itf8(reader)?;

    let number_of_records = read_itf8_as(reader)?;

    let record_counter = read_ltf8(reader).and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let bases = read_ltf8(reader).and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let number_of_blocks = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let landmarks = read_landmarks(reader)?;

    let actual_crc32 = reader.crc().sum();
    let expected_crc32 = reader.get_mut().read_u32::<LittleEndian>()?;

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
        number_of_blocks,
        expected_crc32,
    ) {
        return Ok(0);
    }

    let reference_sequence_context =
        build_reference_sequence_context(reference_sequence_id, alignment_start, alignment_span)?;

    *header = Header::builder()
        .set_reference_sequence_context(reference_sequence_context)
        .set_record_count(number_of_records)
        .set_record_counter(record_counter)
        .set_base_count(bases)
        .set_block_count(number_of_blocks)
        .set_landmarks(landmarks)
        .build();

    Ok(len)
}

fn read_landmarks<R>(reader: &mut R) -> io::Result<Vec<usize>>
where
    R: Read,
{
    let len = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = Vec::with_capacity(len);

    for _ in 0..len {
        let pos = read_itf8(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        buf.push(pos);
    }

    Ok(buf)
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

pub(crate) fn build_reference_sequence_context(
    raw_reference_sequence_id: i32,
    raw_alignment_start: i32,
    raw_alignment_span: i32,
) -> io::Result<ReferenceSequenceContext> {
    const UNMAPPED: i32 = -1;
    const MULTIREF: i32 = -2;

    match raw_reference_sequence_id {
        UNMAPPED => Ok(ReferenceSequenceContext::None),
        MULTIREF => Ok(ReferenceSequenceContext::Many),
        n => {
            let reference_sequence_id =
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let alignment_start = usize::try_from(raw_alignment_start)
                .and_then(Position::try_from)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let alignment_span = usize::try_from(raw_alignment_span)
                .and_then(NonZeroUsize::try_from)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let alignment_end = alignment_start
                .checked_add(usize::from(alignment_span) - 1)
                .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

            Ok(ReferenceSequenceContext::some(
                reference_sequence_id,
                alignment_start,
                alignment_end,
            ))
        }
    }
}

#[cfg(test)]
mod tests {
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

        let expected = Header::builder()
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
