use std::io::{self, Write};

use noodles_vcf::header::string_maps::StringStringMap;

use crate::record::codec::encoder::string_map::write_string_map_index;

pub(super) fn write_key<W>(
    writer: &mut W,
    string_map: &StringStringMap,
    key: &str,
) -> io::Result<()>
where
    W: Write,
{
    string_map
        .get_index_of(key.as_ref())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("genotype key not in string map: {key:?}"),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}
