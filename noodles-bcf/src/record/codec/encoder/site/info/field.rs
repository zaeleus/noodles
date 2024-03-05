mod key;
mod value;

use std::io::{self, Write};

use noodles_vcf::{header::string_maps::StringStringMap, variant::record::info::field::Value};

use self::{key::write_key, value::write_value};

pub(super) fn write_field<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    key: &str,
    value: Option<Value<'_>>,
) -> io::Result<()>
where
    W: Write,
{
    write_key(writer, string_string_map, key)?;
    write_value(writer, value)?;
    Ok(())
}
