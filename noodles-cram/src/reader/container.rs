mod block;
pub(crate) mod header;

pub use self::{block::read_block, header::read_header};

use std::io::{self, Read};

use bytes::BytesMut;

use crate::Container;

pub fn read_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<Option<Container>>
where
    R: Read,
{
    let header = match read_header(reader)? {
        Some(header) => header,
        None => return Ok(None),
    };

    buf.resize(header.len(), 0);
    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(&mut buf)?;
        blocks.push(block);
    }

    Ok(Some(Container::new(header, blocks)))
}
