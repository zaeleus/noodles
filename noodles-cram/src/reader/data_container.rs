pub mod compression_header;
pub mod slice;

pub use self::{compression_header::get_compression_header, slice::read_slice};

use std::io::{self, Read};

use bytes::{Bytes, BytesMut};

use super::container;
use crate::{data_container::CompressionHeader, DataContainer};

pub fn read_data_container<R>(reader: &mut R) -> io::Result<Option<DataContainer>>
where
    R: Read,
{
    let header = container::read_header(reader)?;

    if header.is_eof() {
        return Ok(None);
    }

    let mut buf = BytesMut::new();
    buf.resize(header.len(), 0);
    reader.read_exact(&mut buf)?;
    let mut buf = buf.freeze();

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
) -> io::Result<Option<(crate::container::Header, DataContainer)>>
where
    R: Read,
{
    let header = container::read_header(reader)?;

    if header.is_eof() {
        return Ok(None);
    }

    let mut buf = BytesMut::new();
    buf.resize(header.len(), 0);
    reader.read_exact(&mut buf)?;
    let mut buf = buf.freeze();

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
    let data = block.decompressed_data()?;
    let mut buf = Bytes::from(data.to_vec());
    get_compression_header(&mut buf)
}
