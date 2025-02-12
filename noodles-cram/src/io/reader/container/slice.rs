mod header;

use std::io;

use bytes::Bytes;

use self::header::get_header;
use super::{read_block, Block};
use crate::container::{block::ContentType, Slice};

pub fn read_slice(src: &mut Bytes) -> io::Result<Slice> {
    let header = get_header(src)?;

    let core_data_block = read_core_data_block(src)?;

    let external_block_count = header.block_count() - 1;
    let external_blocks = read_external_blocks(src, external_block_count)?;

    Ok(Slice::new(header, core_data_block, external_blocks))
}

fn read_core_data_block(src: &mut Bytes) -> io::Result<Block> {
    let block = read_block(src)?;

    if block.content_type != ContentType::CoreData {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                ContentType::CoreData,
                block.content_type
            ),
        ));
    }

    Ok(block)
}

fn read_external_blocks(src: &mut Bytes, len: usize) -> io::Result<Vec<Block>> {
    let mut external_blocks = Vec::with_capacity(len);

    for _ in 0..len {
        let block = read_block(src)?;

        if block.content_type != ContentType::ExternalData {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid block content type: expected {:?}, got {:?}",
                    ContentType::ExternalData,
                    block.content_type
                ),
            ));
        }

        external_blocks.push(block);
    }

    Ok(external_blocks)
}
