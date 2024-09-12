mod array;
mod character;
mod genotype;
mod string;

use std::io::{self, Write};

use self::{
    array::write_array, character::write_character, genotype::write_genotype, string::write_string,
};
use crate::{variant::record::samples::series::Value, Header};

pub(super) fn write_value<W>(writer: &mut W, header: &Header, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write_character(writer, *c),
        Value::String(s) => write_string(writer, s),
        Value::Genotype(genotype) => write_genotype(writer, header, genotype.as_ref()),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::*;
    use crate::variant::record_buf::samples::sample::Value as ValueBuf;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, header: &Header, value: Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, header, &value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, Value::Integer(8), b"8")?;
        t(&mut buf, &header, Value::Float(0.333), b"0.333")?;
        t(&mut buf, &header, Value::Character('n'), b"n")?;
        t(
            &mut buf,
            &header,
            Value::String(Cow::from("noodles")),
            b"noodles",
        )?;

        let value_buf = ValueBuf::from(vec![Some(8)]);
        t(&mut buf, &header, (&value_buf).into(), b"8")?;

        Ok(())
    }
}
