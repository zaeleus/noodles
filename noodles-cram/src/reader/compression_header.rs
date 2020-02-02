mod data_series_encodings;
mod preservation_map;
mod tag_encodings;

use std::io::{self, Read};

use self::{
    data_series_encodings::read_data_series_encodings, preservation_map::read_preservation_map,
    tag_encodings::read_tag_encodings,
};

pub fn read_compression_header<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    read_preservation_map(reader)?;
    read_data_series_encodings(reader)?;
    read_tag_encodings(reader)?;
    Ok(())
}
