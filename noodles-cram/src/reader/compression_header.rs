mod data_series_encoding_map;
mod preservation_map;
mod tag_encoding_map;

use std::io::{self, Read};

use self::{
    data_series_encoding_map::read_data_series_encoding_map,
    preservation_map::read_preservation_map, tag_encoding_map::read_tag_encoding_map,
};

use crate::container::CompressionHeader;

pub fn read_compression_header<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: Read,
{
    let preservation_map = read_preservation_map(reader)?;
    let data_series_encoding_map = read_data_series_encoding_map(reader)?;
    let tag_encoding_map = read_tag_encoding_map(reader)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encoding_map,
        tag_encoding_map,
    ))
}
