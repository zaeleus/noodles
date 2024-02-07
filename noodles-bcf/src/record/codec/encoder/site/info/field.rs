mod key;
mod value;

use std::io::{self, Write};

use noodles_vcf as vcf;

use self::{key::write_key, value::write_value};
use crate::header::string_maps::StringStringMap;

pub(super) fn write_field<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    key: &str,
    value: Option<&vcf::variant::record_buf::info::field::Value>,
) -> io::Result<()>
where
    W: Write,
{
    write_key(writer, string_string_map, key)?;
    write_value(writer, value)?;
    Ok(())
}
