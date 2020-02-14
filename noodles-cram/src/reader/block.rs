use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::{num::read_itf8, Block};

pub fn read_block<R>(reader: &mut R, block: &mut Block) -> io::Result<()>
where
    R: Read,
{
    let method = reader.read_u8()?;
    let block_content_type_id = reader.read_u8()?;
    let block_content_id = read_itf8(reader)?;
    let size_in_bytes = read_itf8(reader)?;
    let raw_size_in_bytes = read_itf8(reader)?;

    *block.compression_method_mut() = method;
    *block.content_type_mut() = block_content_type_id;
    *block.content_id_mut() = block_content_id;
    *block.uncompressed_len_mut() = raw_size_in_bytes;

    let data = block.data_mut();
    data.resize(size_in_bytes as usize, Default::default());
    reader.read_exact(data)?;

    let _crc32 = read_crc32(reader)?;

    Ok(())
}

fn read_crc32<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
