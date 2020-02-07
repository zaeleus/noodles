use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::{
    num::{read_itf8, Itf8},
    Block,
};

pub fn read_block<R>(reader: &mut R, block: &mut Block) -> io::Result<()>
where
    R: Read,
{
    let method = read_method(reader)?;
    let content_type_id = read_block_content_type_id(reader)?;
    let content_id = read_block_content_id(reader)?;
    let compressed_len = read_size_in_bytes(reader)?;
    let uncompressed_len = read_raw_size_in_bytes(reader)?;

    *block.compression_method_mut() = method;
    *block.content_type_mut() = content_type_id;
    *block.content_id_mut() = content_id;
    *block.uncompressed_len_mut() = uncompressed_len;

    let data = block.data_mut();
    data.resize(compressed_len as usize, Default::default());
    reader.read_exact(data)?;

    let _crc32 = read_crc32(reader)?;

    Ok(())
}

fn read_method<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    reader.read_u8()
}

fn read_block_content_type_id<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    reader.read_u8()
}

fn read_block_content_id<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_size_in_bytes<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_raw_size_in_bytes<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_crc32<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
