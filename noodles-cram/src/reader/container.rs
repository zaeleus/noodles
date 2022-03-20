mod block;
mod header;

pub use self::{block::read_block, header::read_header};

use std::io::{self, Read};

use bytes::BytesMut;

use crate::Container;

pub fn read_container<R>(reader: &mut R) -> io::Result<Container>
where
    R: Read,
{
    let header = read_header(reader)?;

    let mut buf = BytesMut::new();
    buf.resize(header.len(), 0);
    reader.read_exact(&mut buf)?;
    let mut buf = buf.freeze();

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(&mut buf)?;
        blocks.push(block);
    }

    Ok(Container::new(header, blocks))
}
