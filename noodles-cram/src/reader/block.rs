use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{
    container::{
        block::{CompressionMethod, ContentType},
        Block,
    },
    num::read_itf8,
};

pub fn read_block<R>(reader: &mut R) -> io::Result<Block>
where
    R: Read,
{
    let method = reader.read_u8().and_then(|b| {
        CompressionMethod::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_type_id = reader.read_u8().and_then(|b| {
        ContentType::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_id = read_itf8(reader)?;
    let size_in_bytes = read_itf8(reader)?;
    let raw_size_in_bytes = read_itf8(reader)?;

    let mut data = vec![0; size_in_bytes as usize];
    reader.read_exact(&mut data)?;

    let crc32 = reader.read_u32::<LittleEndian>()?;

    Ok(Block::builder()
        .set_compression_method(method)
        .set_content_type(block_content_type_id)
        .set_content_id(block_content_id)
        .set_uncompressed_len(raw_size_in_bytes)
        .set_data(data)
        .set_crc32(crc32)
        .build())
}
