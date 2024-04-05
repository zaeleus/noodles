mod array;

use std::io::{self, Write};

use self::array::write_array;
use crate::variant::record::info::field::Value;

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Flag => Ok(()),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => writer.write_all(s.as_bytes()),
        Value::Array(array) => write_array(writer, array),
    }
}
