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

use crate::data_container::CompressionHeader;

pub fn get_compression_header(src: &mut Bytes) -> io::Result<CompressionHeader> {
    let preservation_map = get_preservation_map(src)?;
    let data_series_encodings = get_data_series_encodings(src)?;
    let tag_encodings = get_tag_encodings(src)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encodings,
        tag_encodings,
    ))
}
