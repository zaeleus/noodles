mod header;

use bytes::BytesMut;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::header::read_header;
use crate::{
    data_container::DataContainer,
    io::reader::data_container::{read_compression_header_from_block, read_slice},
};

pub async fn read_data_container<R>(
    reader: &mut R,
    buf: &mut BytesMut,
) -> io::Result<Option<DataContainer>>
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

    let compression_header = read_compression_header_from_block(&mut buf)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(&mut buf)?;
        slices.push(slice);
    }

    Ok(Some(DataContainer::new(compression_header, slices)))
}
