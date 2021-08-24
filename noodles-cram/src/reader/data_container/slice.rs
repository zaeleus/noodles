mod header;

pub use self::header::read_header;

use std::io::{self, Read};

use crate::{
    container::Block,
    data_container::{slice, Slice},
    reader::container::read_block,
};

pub fn read_slice<R>(reader: &mut R) -> io::Result<Slice>
where
    R: Read,
{
    let header = read_header_from_block(reader)?;
    let core_data_block = read_block(reader)?;

    let external_block_count = header.block_count() - 1;
    let external_blocks = read_external_blocks(reader, external_block_count)?;

    Ok(Slice::new(header, core_data_block, external_blocks))
}

fn read_header_from_block<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: Read,
{
    let block = read_block(reader)?;
    let data = block.decompressed_data()?;
    let mut data_reader = &data[..];
    read_header(&mut data_reader)
}

fn read_external_blocks<R>(reader: &mut R, len: usize) -> io::Result<Vec<Block>>
where
    R: Read,
{
    let mut external_blocks = Vec::with_capacity(len);

    for _ in 0..len {
        let block = read_block(reader)?;
        external_blocks.push(block);
    }

    Ok(external_blocks)
}
