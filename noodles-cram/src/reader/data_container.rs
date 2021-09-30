pub mod compression_header;
pub mod slice;

pub use self::{compression_header::read_compression_header, slice::read_slice};
use crate::data_container::CompressionHeader;

use std::io::{self, Read};

use super::container;
use crate::DataContainer;

pub fn read_data_container<R>(reader: &mut R) -> io::Result<Option<DataContainer>>
where
    R: Read,
{
    let header = container::read_header(reader)?;

    if header.is_eof() {
        return Ok(None);
    }

    let compression_header = read_compression_header_from_block(reader)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(reader)?;
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

    let compression_header = read_compression_header_from_block(reader)?;

    let slice_count = header.landmarks().len();
    let mut slices = Vec::with_capacity(slice_count);

    for _ in 0..slice_count {
        let slice = read_slice(reader)?;
        slices.push(slice);
    }

    let data_container = DataContainer::new(compression_header, slices);

    Ok(Some((header, data_container)))
}

fn read_compression_header_from_block<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: Read,
{
    use super::container::read_block;

    let block = read_block(reader)?;
    let data = block.decompressed_data()?;
    let mut data_reader = &data[..];
    read_compression_header(&mut data_reader)
}
