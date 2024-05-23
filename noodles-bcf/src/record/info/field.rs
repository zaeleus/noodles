mod value;

use std::io;

use noodles_vcf::{self as vcf, variant::record::info::field::Value};

use self::value::read_value;
use crate::record::codec::decoder::read_string_map_entry;

pub(super) fn read_field<'a, 'h: 'a>(
    src: &mut &'a [u8],
    header: &'h vcf::Header,
) -> io::Result<(&'a str, Option<Value<'a>>)> {
    let key = read_string_map_entry(src, header.string_maps().strings())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let (number, ty) = header
        .infos()
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing info map entry"))?;

    let value = read_value(src, number, ty)?;

    Ok((key, value))
}
