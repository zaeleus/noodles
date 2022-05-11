use noodles_core::Position;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{
    container::{Header, ReferenceSequenceId},
    reader::num::{read_itf8, read_ltf8},
};

// § 9 "End of file container" (2022-04-12)
const EOF_LENGTH: usize = 15;
const EOF_ALIGNMENT_START: usize = 4_542_278;
const EOF_BLOCK_COUNT: usize = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

pub fn read_header<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: Read,
{
    let length = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence_id = read_itf8(reader).and_then(|n| {
        ReferenceSequenceId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let starting_position_on_the_reference = read_itf8(reader)
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(Position::new)?;

    let alignment_span = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let number_of_records = read_itf8(reader)?;
    let record_counter = read_ltf8(reader)?;

    let bases = read_ltf8(reader).and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let number_of_blocks = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let landmarks = read_landmarks(reader)?;
    let crc32 = reader.read_u32::<LittleEndian>()?;

    if is_eof(
        length,
        reference_sequence_id,
        starting_position_on_the_reference,
        number_of_blocks,
        crc32,
    ) {
        return Ok(None);
    }

    let mut builder = Header::builder()
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

    Ok(Some(builder.build()))
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
    reference_sequence_id: ReferenceSequenceId,
    alignment_start: Option<Position>,
    block_count: usize,
    crc32: u32,
) -> bool {
    length == EOF_LENGTH
        && reference_sequence_id.is_none()
        && alignment_start
            .map(|position| usize::from(position) == EOF_ALIGNMENT_START)
            .unwrap_or(false)
        && block_count == EOF_BLOCK_COUNT
        && crc32 == EOF_CRC32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
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
        let actual = read_header(&mut reader)?;

        let expected = Header::builder()
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

        assert_eq!(actual, Some(expected));

        Ok(())
    }

    #[test]
    fn test_read_header_with_eof() -> io::Result<()> {
        let data = [
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
        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        assert!(actual.is_none());

        Ok(())
    }
}
