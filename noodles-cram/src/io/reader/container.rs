mod block;
pub mod compression_header;
pub mod header;
pub mod slice;

use std::io::{self, Read};

use bytes::BytesMut;

use self::header::read_header;
pub use self::{block::read_block, compression_header::get_compression_header, slice::read_slice};
use crate::{container::Header, Container};

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

    let compression_header = get_compression_header(&mut buf)?;

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

    let compression_header = get_compression_header(&mut buf)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(&mut buf)?;
        slices.push(slice);
    }

    let container = Container::new(compression_header, slices);

    Ok(Some((header, len, container)))
}
