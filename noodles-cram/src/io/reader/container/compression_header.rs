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
use super::read_block;
use crate::container::{block::ContentType, CompressionHeader};

pub fn read_compression_header(src: &mut &[u8]) -> io::Result<CompressionHeader> {
    let block = read_block(src)?;

    if block.content_type != ContentType::CompressionHeader {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                ContentType::CompressionHeader,
                block.content_type
            ),
        ));
    }

    let buf = block.decode()?;
    read_compression_header_inner(&mut &buf[..])
}

fn read_compression_header_inner(src: &mut &[u8]) -> io::Result<CompressionHeader> {
    let preservation_map = read_preservation_map(src)?;
    let data_series_encodings = read_data_series_encodings(src)?;
    let tag_encodings = read_tag_encodings(src)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encodings,
        tag_encodings,
    ))
}
