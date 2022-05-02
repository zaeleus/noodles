use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;
use noodles_core::Position;

use crate::{
    container,
    writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &container::Header) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length =
        i32::try_from(header.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    crc_writer.write_i32::<LittleEndian>(length)?;

    let reference_sequence_id = i32::try_from(header.reference_sequence_id())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, reference_sequence_id)?;

    write_starting_position_on_the_reference(&mut crc_writer, header.start_position())?;

    let alignment_span = i32::try_from(header.alignment_span())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, alignment_span)?;

    let number_of_records = header.record_count();
    write_itf8(&mut crc_writer, number_of_records)?;

    let record_counter = header.record_counter();
    write_ltf8(&mut crc_writer, record_counter)?;

    let bases = i64::try_from(header.base_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(&mut crc_writer, bases)?;

    let number_of_blocks = i32::try_from(header.block_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, number_of_blocks)?;

    write_landmarks(&mut crc_writer, header.landmarks())?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}

fn write_starting_position_on_the_reference<W>(
    writer: &mut W,
    start_position: Option<Position>,
) -> io::Result<()>
where
    W: Write,
{
    let n = start_position.map(usize::from).unwrap_or_default();

    let starting_position_on_the_reference =
        i32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_itf8(writer, starting_position_on_the_reference)
}

fn write_landmarks<W>(writer: &mut W, landmarks: &[usize]) -> io::Result<()>
where
    W: Write,
{
    let landmarks_len = i32::try_from(landmarks.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, landmarks_len)?;

    for &pos in landmarks {
        let n = i32::try_from(pos).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_itf8(writer, n)?;
    }

    Ok(())
}
