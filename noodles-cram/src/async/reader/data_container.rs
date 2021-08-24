pub mod compression_header;
pub mod slice;

pub use self::{compression_header::read_compression_header, slice::read_slice};

use tokio::io::{self, AsyncRead};

use crate::data_container::{CompressionHeader, DataContainer};

pub async fn read_data_container<R>(reader: &mut R) -> io::Result<Option<DataContainer>>
where
    R: AsyncRead + Unpin,
{
    let header = super::container::read_header(reader).await?;

    if header.is_eof() {
        return Ok(None);
    }

    let compression_header = read_compression_header_from_block(reader).await?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(reader).await?;
        slices.push(slice);
    }

    Ok(Some(DataContainer::new(compression_header, slices)))
}

async fn read_compression_header_from_block<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: AsyncRead + Unpin,
{
    use super::container::read_block;

    let block = read_block(reader).await?;
    let data = block.decompressed_data()?;
    let mut data_reader = &data[..];
    read_compression_header(&mut data_reader).await
}
