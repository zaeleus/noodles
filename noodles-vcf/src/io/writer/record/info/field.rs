mod value;

use std::io::{self, Write};

use self::value::write_value;
use crate::{io::writer::record::MISSING, variant::record::info::field::Value};

pub(super) fn write_field<W>(writer: &mut W, key: &str, value: Option<&Value>) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b"=";

    writer.write_all(key.as_bytes())?;

    match value {
        Some(Value::Flag) => {}
        Some(v) => {
            writer.write_all(SEPARATOR)?;
            write_value(writer, v)?;
        }
        None => {
            writer.write_all(SEPARATOR)?;
            writer.write_all(MISSING)?;
        }
    }

    Ok(())
}
