use std::io;

use noodles_vcf as vcf;

use crate::{
    header::string_maps::StringStringMap, record::codec::decoder::string_map::read_string_map_index,
};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    infos: &'h vcf::header::Infos,
    string_string_map: &StringStringMap,
) -> io::Result<&'h vcf::record::info::field::Key> {
    read_string_map_index(src)
        .and_then(|j| {
            string_string_map.get_index(j).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid string map index: {j}"),
                )
            })
        })
        .and_then(|raw_key| {
            infos.get_key_value(raw_key).map(|(k, _)| k).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("missing header INFO record for {raw_key}"),
                )
            })
        })
}
