mod data_series_encodings;
mod encoding;
mod preservation_map;
mod tag_encodings;

use std::io;

use bytes::Bytes;

use self::{
    data_series_encodings::get_data_series_encodings,
    encoding::{
        get_encoding_for_byte_array_codec, get_encoding_for_byte_codec,
        get_encoding_for_integer_codec,
    },
    preservation_map::get_preservation_map,
    tag_encodings::get_tag_encodings,
};
use super::read_block;
use crate::container::{block::ContentType, CompressionHeader};

pub fn get_compression_header(src: &mut Bytes) -> io::Result<CompressionHeader> {
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

    let mut buf = block.decode()?;
    get_compression_header_inner(&mut buf)
}

fn get_compression_header_inner(src: &mut Bytes) -> io::Result<CompressionHeader> {
    let preservation_map = get_preservation_map(src)?;
    let data_series_encodings = get_data_series_encodings(src)?;
    let tag_encodings = get_tag_encodings(src)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encodings,
        tag_encodings,
    ))
}
