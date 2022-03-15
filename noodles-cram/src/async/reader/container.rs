mod header;

pub use self::header::read_header;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{reader::container::read_block, Container};

pub async fn read_container<R>(reader: &mut R) -> io::Result<Container>
where
    R: AsyncRead + Unpin,
{
    let header = read_header(reader).await?;

    let mut buf = vec![0; header.len()];
    reader.read_exact(&mut buf).await?;

    let mut buf_reader = &buf[..];

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(&mut buf_reader)?;
        blocks.push(block);
    }

    Ok(Container::new(header, blocks))
}
