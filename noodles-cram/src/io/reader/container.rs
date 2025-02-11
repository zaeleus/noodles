mod block;
pub mod compression_header;
pub mod header;
pub mod slice;

use std::io::{self, Read};

use bytes::{Bytes, BytesMut};

use self::header::read_header;
pub use self::{block::read_block, compression_header::get_compression_header, slice::read_slice};
use crate::{
    container::{block::ContentType, CompressionHeader, Header},
    Container,
};

pub fn read_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<Option<Container>>
where
    R: Read,
{
    let mut header = Header::default();

    let len = match read_header(reader, &mut header)? {
        0 => return Ok(None),
        n => n,
    };

    buf.resize(len, 0);
    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    let compression_header = read_compression_header_from_block(&mut buf)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(&mut buf)?;
        slices.push(slice);
    }

    Ok(Some(Container::new(compression_header, slices)))
}

pub fn read_container_with_header<R>(
    reader: &mut R,
    buf: &mut BytesMut,
) -> io::Result<Option<(crate::container::Header, usize, Container)>>
where
    R: Read,
{
    let mut header = Header::default();

    let len = match read_header(reader, &mut header)? {
        0 => return Ok(None),
        n => n,
    };

    buf.resize(len, 0);

    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    let compression_header = read_compression_header_from_block(&mut buf)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(&mut buf)?;
        slices.push(slice);
    }

    let container = Container::new(compression_header, slices);

    Ok(Some((header, len, container)))
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
