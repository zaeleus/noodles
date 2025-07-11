use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use flate2::CrcReader;

use crate::io::reader::num::{read_i32_le, read_itf8, read_itf8_as, read_ltf8};

pub(crate) fn read_header<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut crc_reader = CrcReader::new(reader);
    read_header_inner(&mut crc_reader)
}

fn read_header_inner<R>(reader: &mut CrcReader<R>) -> io::Result<u64>
where
    R: Read,
{
    let length = read_i32_le(reader).and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let _reference_sequence_id = read_itf8(reader)?;
    let _alignment_start = read_itf8(reader)?;
    let _alignment_span = read_itf8(reader)?;
    let _record_count = read_itf8(reader)?;
    let _record_counter = read_ltf8(reader)?;
    let _base_count = read_ltf8(reader)?;
    let _block_count = read_itf8(reader)?;
    read_landmarks(reader)?;

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

    Ok(length)
}

fn read_landmarks<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let len: usize = read_itf8_as(reader)?;

    for _ in 0..len {
        read_itf8(reader)?;
    }

    Ok(())
}
