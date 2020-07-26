mod data_series_encoding_map;
mod preservation_map;
mod tag_encoding_map;

use std::io::{self, Read};

use self::{
    data_series_encoding_map::read_data_series_encoding_map,
    preservation_map::read_preservation_map, tag_encoding_map::read_tag_encoding_map,
};

use crate::container::CompressionHeader;

use super::block::read_block;

pub fn read_compression_header<R>(reader: &mut R) -> io::Result<CompressionHeader>
where
    R: Read,
{
    let block = read_block(reader)?;

    let data = block.decompressed_data();
    let mut data_reader = &data[..];

    let preservation_map = read_preservation_map(&mut data_reader)?;
    let data_series_encoding_map = read_data_series_encoding_map(&mut data_reader)?;
    let tag_encoding_map = read_tag_encoding_map(&mut data_reader)?;

    Ok(CompressionHeader::new(
        preservation_map,
        data_series_encoding_map,
        tag_encoding_map,
    ))
}
