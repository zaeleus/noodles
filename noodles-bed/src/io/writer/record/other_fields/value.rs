mod character;

use std::io::{self, Write};

use self::character::write_character;
use crate::feature::record::other_fields::Value;

pub(super) fn write_value<W>(writer: &mut W, value: Value<'_>) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Int64(n) => write!(writer, "{n}"),
        Value::UInt64(n) => write!(writer, "{n}"),
        Value::Float64(n) => write!(writer, "{n}"),
        Value::Character(b) => write_character(writer, b),
        Value::String(s) => writer.write_all(s),
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: Value<'_>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Value::Int64(-8), b"-8")?;
        t(&mut buf, Value::UInt64(13), b"13")?;
        t(&mut buf, Value::Float64(0.0), b"0")?;
        t(&mut buf, Value::Character(b'n'), b"n")?;
        t(&mut buf, Value::String(b"ndls".as_bstr()), b"ndls")?;

        Ok(())
    }
}
