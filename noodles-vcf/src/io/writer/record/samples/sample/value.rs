mod array;
mod character;
mod genotype;
mod integer;
mod string;

use std::io::{self, Write};

use self::{
    array::write_array, character::write_character, genotype::write_genotype,
    integer::write_integer, string::write_string,
};
use crate::{Header, variant::record::samples::series::Value};

pub(super) fn write_value<W>(writer: &mut W, header: &Header, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write_integer(writer, *n),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write_character(writer, *c),
        Value::String(s) => write_string(writer, s),
        Value::Genotype(genotype) => write_genotype(writer, header, genotype.as_ref()),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::samples::sample::Value as ValueBuf;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            value: &ValueBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_value(buf, header, &Value::from(value))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, &ValueBuf::from(8), b"8")?;
        t(&mut buf, &header, &ValueBuf::from(0.333), b"0.333")?;
        t(&mut buf, &header, &ValueBuf::from('n'), b"n")?;
        t(&mut buf, &header, &ValueBuf::from("noodles"), b"noodles")?;
        t(&mut buf, &header, &ValueBuf::from(vec![Some(8)]), b"8")?;

        Ok(())
    }
}
