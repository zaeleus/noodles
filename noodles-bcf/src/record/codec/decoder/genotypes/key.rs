use std::io;

use noodles_vcf::{self as vcf, record::genotypes::keys::Key};

use crate::{
    header::string_maps::StringStringMap, record::codec::decoder::string_map::read_string_map_index,
};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    formats: &'h vcf::header::Formats,
    string_map: &StringStringMap,
) -> io::Result<&'h Key> {
    read_string_map_index(src)
        .and_then(|j| {
            string_map.get_index(j).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid string map index: {j}"),
                )
            })
        })
        .and_then(|raw_key| {
            formats
                .get_key_value(raw_key)
                .map(|(k, _)| k)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing header FORMAT record for {raw_key}"),
                    )
                })
        })
}
