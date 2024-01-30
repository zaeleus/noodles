use std::io::{self, Write};

use crate::{
    header::string_maps::StringStringMap,
    record::codec::encoder::string_map::write_string_map_index,
};

pub(super) fn write_key<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    key: &str,
) -> io::Result<()>
where
    W: Write,
{
    string_string_map
        .get_index_of(key)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("info key missing from string map: {key:?}"),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}
