mod header;

pub use self::header::get_header;

use std::io;

use bytes::Bytes;

use crate::{
    container::Block,
    data_container::{slice, Slice},
    reader::container::read_block,
};

pub fn read_slice(src: &mut Bytes) -> io::Result<Slice> {
    let header = read_header_from_block(src)?;
    let core_data_block = read_block(src)?;

    let external_block_count = header.block_count() - 1;
    let external_blocks = read_external_blocks(src, external_block_count)?;

    Ok(Slice::new(header, core_data_block, external_blocks))
}

fn read_header_from_block(src: &mut Bytes) -> io::Result<slice::Header> {
    let block = read_block(src)?;
    let mut data = block.decompressed_data()?;
    get_header(&mut data)
}

fn read_external_blocks(src: &mut Bytes, len: usize) -> io::Result<Vec<Block>> {
    let mut external_blocks = Vec::with_capacity(len);

    for _ in 0..len {
        let block = read_block(src)?;
        external_blocks.push(block);
    }

    Ok(external_blocks)
}
