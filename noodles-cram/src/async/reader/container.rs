mod header;

pub use self::header::read_header;

use bytes::BytesMut;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{reader::container::read_block, Container};

pub async fn read_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<Option<Container>>
where
    R: AsyncRead + Unpin,
{
    let header = match read_header(reader).await? {
        Some(header) => header,
        None => return Ok(None),
    };

    buf.resize(header.len(), 0);
    reader.read_exact(buf).await?;
    let mut buf = buf.split().freeze();

    let blocks_len = header.block_count();
    let mut blocks = Vec::with_capacity(blocks_len);

    for _ in 0..blocks_len {
        let block = read_block(&mut buf)?;
        blocks.push(block);
    }

    Ok(Some(Container::new(header, blocks)))
}
