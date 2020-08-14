mod data_series_encoding_map;
mod preservation_map;
mod tag_encoding_map;

use std::io::{self, Write};

use self::{
    data_series_encoding_map::write_data_series_encoding_map,
    preservation_map::write_preservation_map, tag_encoding_map::write_tag_encoding_map,
};

use crate::container::CompressionHeader;

#[allow(dead_code)]
pub fn write_compression_header<W>(
    writer: &mut W,
    compression_header: &CompressionHeader,
) -> io::Result<()>
where
    W: Write,
{
    write_preservation_map(writer, compression_header.preservation_map())?;
    write_data_series_encoding_map(writer, compression_header.data_series_encoding_map())?;
    write_tag_encoding_map(writer, compression_header.tag_encoding_map())?;
    Ok(())
}
