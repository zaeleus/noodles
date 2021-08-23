mod header;

pub use self::header::read_header;

use tokio::io::{self, AsyncRead};

use super::container::read_block;
use crate::{
    container::Block,
    data_container::{slice, Slice},
};

pub async fn read_slice<R>(reader: &mut R) -> io::Result<Slice>
where
    R: AsyncRead + Unpin,
{
    let header = read_header_from_block(reader).await?;
    let core_data_block = read_block(reader).await?;

    let external_block_count = header.block_count() - 1;
    let external_blocks = read_external_blocks(reader, external_block_count).await?;

    Ok(Slice::new(header, core_data_block, external_blocks))
}

async fn read_header_from_block<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: AsyncRead + Unpin,
{
    let header_block = read_block(reader).await?;
    let data = header_block.decompressed_data()?;
    let mut data_reader = &data[..];
    read_header(&mut data_reader).await
}

async fn read_external_blocks<R>(reader: &mut R, len: usize) -> io::Result<Vec<Block>>
where
    R: AsyncRead + Unpin,
{
    let mut external_blocks = Vec::with_capacity(len);

    for _ in 0..len {
        let block = read_block(reader).await?;
        external_blocks.push(block);
    }

    Ok(external_blocks)
}
