mod data_series_encoding_map;
mod encoding;
mod preservation_map;
mod tag_encoding_map;

use std::io;

use bytes::Bytes;

use self::{
    data_series_encoding_map::get_data_series_encoding_map,
    encoding::{
        get_encoding_for_byte_array_codec, get_encoding_for_byte_codec,
        get_encoding_for_integer_codec,
    },
    preservation_map::get_preservation_map,
    tag_encoding_map::get_tag_encoding_map,
};

use crate::data_container::CompressionHeader;

pub fn get_compression_header(src: &mut Bytes) -> io::Result<CompressionHeader> {
    let preservation_map = get_preservation_map(src)?;
    let data_series_encoding_map = get_data_series_encoding_map(src)?;
    let tag_encoding_map = get_tag_encoding_map(src)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encoding_map,
        tag_encoding_map,
    ))
}
