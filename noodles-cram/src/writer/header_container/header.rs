use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use crate::writer::num::{write_itf8, write_ltf8};

pub(super) fn write_header<W>(writer: &mut W, len: usize) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    crc_writer.write_i32::<LittleEndian>(length)?;

    // reference sequence ID
    write_itf8(&mut crc_writer, -1)?;

    // alignment start
    write_itf8(&mut crc_writer, 0)?;

    // alignment span
    write_itf8(&mut crc_writer, 0)?;

    // record count
    write_itf8(&mut crc_writer, 0)?;

    // record count
    write_ltf8(&mut crc_writer, 0)?;

    // base count
    write_ltf8(&mut crc_writer, 0)?;

    // block count
    write_itf8(&mut crc_writer, 1)?;

    // landmarks
    write_itf8(&mut crc_writer, 0)?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}
