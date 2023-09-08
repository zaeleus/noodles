mod key;
mod value;

use std::io;

use noodles_vcf as vcf;

use self::{key::read_key, value::read_value};
use crate::header::string_maps::StringStringMap;

pub(crate) fn read_field(
    src: &mut &[u8],
    infos: &vcf::header::Infos,
    string_string_map: &StringStringMap,
) -> io::Result<(
    vcf::record::info::field::Key,
    Option<vcf::record::info::field::Value>,
)> {
    let key = read_key(src, infos, string_string_map)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let info = infos.get(key).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("missing header INFO record for {key}"),
        )
    })?;

    let value = read_value(src, info).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((key.clone(), value))
}
