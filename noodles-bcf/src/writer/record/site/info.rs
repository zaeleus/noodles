use std::io::{self, Write};

use noodles_vcf as vcf;

use crate::{
    header::StringMap,
    record::value::{Float, Int32, Value},
    writer::{string_map::write_string_map_index, value::write_value},
};

pub fn write_info<W>(
    writer: &mut W,
    string_map: &StringMap,
    info: &vcf::record::Info,
) -> io::Result<()>
where
    W: Write,
{
    for field in info.values() {
        write_info_key(writer, string_map, field.key())?;
        write_info_value(writer, field.value())?;
    }

    Ok(())
}

fn write_info_key<W>(
    writer: &mut W,
    string_map: &StringMap,
    key: &vcf::record::info::field::Key,
) -> io::Result<()>
where
    W: Write,
{
    string_map
        .get_index_of(key.as_ref())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("info key missing from string map: {:?}", key),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}

fn write_info_value<W>(writer: &mut W, value: &vcf::record::info::field::Value) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::info::field;

    match value {
        field::Value::Integer(n) => write_value(writer, Some(Value::Int32(Some(Int32::Value(*n))))),
        field::Value::Float(n) => write_value(writer, Some(Value::Float(Some(Float::Value(*n))))),
        field::Value::Flag => write_value(writer, None),
        field::Value::String(s) => write_value(writer, Some(Value::String(Some(s.into())))),
        v => todo!("unhandled INFO field value: {:?}", v),
    }
}
