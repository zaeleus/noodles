mod key;
mod value;

use std::io::{self, Write};

use self::{key::write_key, value::write_value};
use crate::{io::writer::record::MISSING, variant::record::info::field::Value};

pub(super) fn write_field<W>(writer: &mut W, key: &str, value: Option<&Value>) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b"=";

    write_key(writer, key).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    match value {
        Some(Value::Flag) => {}
        Some(v) => {
            writer.write_all(SEPARATOR)?;
            write_value(writer, v).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        }
        None => {
            writer.write_all(SEPARATOR)?;
            writer.write_all(MISSING)?;
        }
    }

    Ok(())
}
