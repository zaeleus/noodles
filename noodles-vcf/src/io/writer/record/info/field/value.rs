mod array;
mod character;
mod float;
mod integer;
mod string;

use std::io::{self, Write};

use self::{
    array::write_array, character::write_character, float::write_float, integer::write_integer,
    string::write_string,
};
use crate::variant::record::info::field::Value;

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write_integer(writer, *n),
        Value::Float(n) => write_float(writer, *n),
        Value::Flag => Ok(()),
        Value::Character(c) => write_character(writer, *c),
        Value::String(s) => write_string(writer, s),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        use crate::variant::record_buf::info::field::Value as ValueBuf;

        fn t(buf: &mut Vec<u8>, value: &ValueBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, &Value::from(value))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &ValueBuf::from(0), b"0")?;
        t(&mut buf, &ValueBuf::from(0.0), b"0")?;
        t(&mut buf, &ValueBuf::Flag, b"")?;
        t(&mut buf, &ValueBuf::from('n'), b"n")?;
        t(&mut buf, &ValueBuf::from("noodles"), b"noodles")?;
        t(&mut buf, &ValueBuf::from(vec![Some(8)]), b"8")?;

        Ok(())
    }
}
