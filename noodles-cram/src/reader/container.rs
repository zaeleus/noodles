mod header;

pub use self::header::read_header;

use std::io::{self, Read};

use crate::Container;

pub fn read_container<R>(reader: &mut R) -> io::Result<Container>
where
    R: Read,
{
    use super::block::read_block;

    let header = read_header(reader)?;

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(reader)?;
        blocks.push(block);
    }

    Ok(Container::new(header, blocks))
}
