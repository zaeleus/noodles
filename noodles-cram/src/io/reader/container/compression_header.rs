mod data_series_encodings;
mod encoding;
mod preservation_map;
mod tag_encodings;

use std::io;

use self::{
    data_series_encodings::read_data_series_encodings,
    encoding::{read_byte_array_encoding, read_byte_encoding, read_integer_encoding},
    preservation_map::read_preservation_map,
    tag_encodings::read_tag_encodings,
};
use super::read_block_as;
use crate::{
    container::{CompressionHeader, block::ContentType},
    file_definition::Version,
};

pub fn read_compression_header(src: &mut &[u8], version: Version) -> io::Result<CompressionHeader> {
    let block = read_block_as(src, ContentType::CompressionHeader, version)?;
    let buf = block.decode()?;
    read_compression_header_inner(&mut &buf[..], version)
}

fn read_compression_header_inner(
    src: &mut &[u8],
    version: Version,
) -> io::Result<CompressionHeader> {
    let preservation_map = read_preservation_map(src, version)?;
    let data_series_encodings = read_data_series_encodings(src, version)?;
    let tag_encodings = read_tag_encodings(src, version)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encodings,
        tag_encodings,
    ))
}
