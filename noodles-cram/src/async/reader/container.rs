mod header;

pub use self::header::read_header;

use bytes::BytesMut;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{reader::container::read_block, Container};

pub async fn read_container<R>(reader: &mut R) -> io::Result<Container>
where
    R: AsyncRead + Unpin,
{
    let header = read_header(reader).await?;

    let mut buf = BytesMut::new();
    buf.resize(header.len(), 0);
    reader.read_exact(&mut buf).await?;
    let mut buf = buf.freeze();

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(&mut buf)?;
        blocks.push(block);
    }

    Ok(Container::new(header, blocks))
}
