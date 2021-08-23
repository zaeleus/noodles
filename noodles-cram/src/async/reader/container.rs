mod block;
mod header;

pub use self::{block::read_block, header::read_header};

use tokio::io::{self, AsyncRead};

use crate::Container;

pub async fn read_container<R>(reader: &mut R) -> io::Result<Container>
where
    R: AsyncRead + Unpin,
{
    let header = read_header(reader).await?;

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(reader).await?;
        blocks.push(block);
    }

    Ok(Container::new(header, blocks))
}
