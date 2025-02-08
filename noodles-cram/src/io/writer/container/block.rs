mod compression_method;
mod content_type;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use self::{compression_method::write_compression_method, content_type::write_content_type};
use crate::{container::Block, io::writer::num::write_itf8};

pub fn write_block<W>(writer: &mut W, block: &Block) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    write_compression_method(&mut crc_writer, block.compression_method())?;
    write_content_type(&mut crc_writer, block.content_type())?;
    write_itf8(&mut crc_writer, block.content_id())?;

    let size_in_bytes = i32::try_from(block.data().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, size_in_bytes)?;

    let uncompressed_data_len = i32::try_from(block.uncompressed_len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, uncompressed_data_len)?;

    crc_writer.write_all(block.data())?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}
