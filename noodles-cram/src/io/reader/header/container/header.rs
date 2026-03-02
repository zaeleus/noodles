use std::io::{self, Read};

use flate2::CrcReader;

use crate::{
    file_definition::Version,
    io::reader::num::{
        read_header_int, read_i32_le, read_long_as, read_position, read_u32_le, read_uint7,
        read_unsigned_int_as,
    },
};

pub(crate) fn read_header<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: Read,
{
    if version.has_crc32() {
        let mut crc_reader = CrcReader::new(reader);
        read_header_with_crc32(&mut crc_reader, version)
    } else {
        read_header_without_crc32(reader, version)
    }
}

fn read_length<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: Read,
{
    if version.uses_vlq() {
        read_uint7(reader).map(u64::from)
    } else {
        read_i32_le(reader).and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}

fn read_header_with_crc32<R>(reader: &mut CrcReader<R>, version: Version) -> io::Result<u64>
where
    R: Read,
{
    let length = read_length(reader, version)?;

    let _reference_sequence_id = read_header_int(reader, version)?;
    let _alignment_start = read_position(reader, version)?;
    let _alignment_span = read_position(reader, version)?;
    let _record_count: usize = read_unsigned_int_as(reader, version)?;
    let _record_counter: usize = read_long_as(reader, version)?;
    let _base_count: usize = read_long_as(reader, version)?;
    let _block_count: usize = read_unsigned_int_as(reader, version)?;
    read_landmarks(reader, version)?;

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

    Ok(length)
}

fn read_header_without_crc32<R>(reader: &mut R, version: Version) -> io::Result<u64>
where
    R: Read,
{
    let length = read_length(reader, version)?;

    let _reference_sequence_id = read_header_int(reader, version)?;
    let _alignment_start = read_position(reader, version)?;
    let _alignment_span = read_position(reader, version)?;
    let _record_count: usize = read_unsigned_int_as(reader, version)?;
    let _record_counter: usize = read_long_as(reader, version)?;
    let _base_count: usize = read_long_as(reader, version)?;
    let _block_count: usize = read_unsigned_int_as(reader, version)?;
    read_landmarks(reader, version)?;

    Ok(length)
}

fn read_landmarks<R>(reader: &mut R, version: Version) -> io::Result<()>
where
    R: Read,
{
    let len: usize = read_unsigned_int_as(reader, version)?;

    for _ in 0..len {
        let _landmark: usize = read_unsigned_int_as(reader, version)?;
    }

    Ok(())
}
