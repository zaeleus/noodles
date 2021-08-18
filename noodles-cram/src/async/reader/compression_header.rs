mod data_series_encoding_map;
mod encoding;
mod preservation_map;
mod tag_encoding_map;

use tokio::io::{self, AsyncRead};

use self::{
    data_series_encoding_map::read_data_series_encoding_map, encoding::read_encoding,
    preservation_map::read_preservation_map, tag_encoding_map::read_tag_encoding_map,
};
use crate::container::CompressionHeader;

pub async fn read_compression_header<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: AsyncRead + Unpin,
{
    let preservation_map = read_preservation_map(reader).await?;
    let data_series_encoding_map = read_data_series_encoding_map(reader).await?;
    let tag_encoding_map = read_tag_encoding_map(reader).await?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encoding_map,
        tag_encoding_map,
    ))
}
