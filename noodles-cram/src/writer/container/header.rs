use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use crate::{
    container,
    writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &container::Header) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length = header.len();
    crc_writer.write_i32::<LittleEndian>(length)?;

    let reference_sequence_id = i32::from(header.reference_sequence_id());
    write_itf8(&mut crc_writer, reference_sequence_id)?;

    let starting_position_on_the_reference =
        header.start_position().map(i32::from).unwrap_or_default();
    write_itf8(&mut crc_writer, starting_position_on_the_reference)?;

    let alignment_span = i32::try_from(header.alignment_span())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, alignment_span)?;

    let number_of_records = header.record_count();
    write_itf8(&mut crc_writer, number_of_records)?;

    let record_counter = header.record_counter();
    write_ltf8(&mut crc_writer, record_counter)?;

    let bases = header.base_count();
    write_ltf8(&mut crc_writer, bases)?;

    let number_of_blocks = i32::try_from(header.block_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, number_of_blocks)?;

    let landmarks_len = header.landmarks().len() as i32;
    write_itf8(&mut crc_writer, landmarks_len)?;

    for &pos in header.landmarks() {
        write_itf8(&mut crc_writer, pos)?;
    }

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}
