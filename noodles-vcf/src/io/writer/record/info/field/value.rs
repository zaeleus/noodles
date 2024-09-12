mod array;
mod character;
mod string;

use std::io::{self, Write};

use self::{array::write_array, character::write_character, string::write_string};
use crate::variant::record::info::field::Value;

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Flag => Ok(()),
        Value::Character(c) => write_character(writer, *c),
        Value::String(s) => write_string(writer, s),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        use crate::variant::record_buf::info::field::Value as ValueBuf;

        fn t(buf: &mut Vec<u8>, value: &Value<'_>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Integer(0), b"0")?;
        t(&mut buf, &Value::Float(0.0), b"0")?;
        t(&mut buf, &Value::Flag, b"")?;

        t(&mut buf, &Value::Character('n'), b"n")?;
        t(&mut buf, &Value::Character(';'), b";")?; // FIXME

        t(&mut buf, &Value::String(Cow::from("ndls")), b"ndls")?;
        t(&mut buf, &Value::String(Cow::from("n;d")), b"n%3Bd")?;

        let value_buf = ValueBuf::from(vec![Some(8), Some(13), None]);
        t(&mut buf, &Value::from(&value_buf), b"8,13,.")?;

        Ok(())
    }
}
