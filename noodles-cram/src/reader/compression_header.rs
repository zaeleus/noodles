mod data_series_encodings;
mod preservation_map;
mod tag_encodings;

use std::io::{self, Read};

use self::{
    data_series_encodings::read_data_series_encodings, preservation_map::read_preservation_map,
    tag_encodings::read_tag_encodings,
};

use crate::CompressionHeader;

pub fn read_compression_header<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: Read,
{
    let preservation_map = read_preservation_map(reader)?;
    let data_series_encodings = read_data_series_encodings(reader)?;
    let tag_encodings = read_tag_encodings(reader)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encodings,
        tag_encodings,
    ))
}
