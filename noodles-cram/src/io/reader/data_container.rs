pub mod compression_header;
pub mod header;
pub mod slice;

pub use self::{compression_header::get_compression_header, slice::read_slice};

use std::io::{self, Read};

use bytes::{Bytes, BytesMut};

use self::header::read_header;
use crate::{container::block::ContentType, data_container::CompressionHeader, DataContainer};

pub fn read_data_container<R>(
    reader: &mut R,
    buf: &mut BytesMut,
) -> io::Result<Option<DataContainer>>
where
    R: Read,
{
    let Some(header) = read_header(reader)? else {
        return Ok(None);
    };

    buf.resize(header.len(), 0);
    reader.read_exact(buf)?;
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

pub fn read_data_container_with_container_header<R>(
    reader: &mut R,
    buf: &mut BytesMut,
) -> io::Result<Option<(crate::data_container::Header, DataContainer)>>
where
    R: Read,
{
    let Some(header) = read_header(reader)? else {
        return Ok(None);
    };

    buf.resize(header.len(), 0);
    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    let compression_header = read_compression_header_from_block(&mut buf)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(&mut buf)?;
        slices.push(slice);
    }

    let data_container = DataContainer::new(compression_header, slices);

    Ok(Some((header, data_container)))
}

pub(crate) fn read_compression_header_from_block(src: &mut Bytes) -> io::Result<CompressionHeader> {
    use super::container::read_block;

    let block = read_block(src)?;

    if block.content_type() != ContentType::CompressionHeader {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                ContentType::CompressionHeader,
                block.content_type()
            ),
        ));
    }

    let mut data = block.decompressed_data()?;
    get_compression_header(&mut data)
}
