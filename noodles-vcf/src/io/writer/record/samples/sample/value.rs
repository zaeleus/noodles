mod array;
mod genotype;
mod string;

use std::io::{self, Write};

use self::{array::write_array, genotype::write_genotype, string::write_string};
use crate::{variant::record::samples::series::Value, Header};

pub(super) fn write_value<W>(writer: &mut W, header: &Header, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => write_string(writer, s),
        Value::Genotype(genotype) => write_genotype(writer, header, genotype.as_ref()),
        Value::Array(array) => write_array(writer, array),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::samples::series::value::Array;

    #[test]
    fn test_write_value() -> io::Result<()> {
        struct Values<'a, T>(&'a [Option<T>]);

        impl<'a, T> crate::variant::record::samples::series::value::array::Values<'a, T> for Values<'a, T>
        where
            T: Copy,
        {
            fn len(&self) -> usize {
                self.0.len()
            }

            fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<T>>> + '_> {
                Box::new(self.0.iter().copied().map(Ok))
            }
        }

        impl<'a> crate::variant::record::samples::series::value::array::Values<'a, &'a str>
            for Values<'a, String>
        {
            fn len(&self) -> usize {
                self.0.len()
            }

            fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<&'a str>>> + '_> {
                Box::new(self.0.iter().map(|s| Ok(s.as_deref())))
            }
        }

        fn t(buf: &mut Vec<u8>, header: &Header, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, header, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, &Value::Integer(8), b"8")?;
        t(&mut buf, &header, &Value::Float(0.333), b"0.333")?;
        t(&mut buf, &header, &Value::Character('n'), b"n")?;
        t(&mut buf, &header, &Value::String("noodles"), b"noodles")?;

        t(
            &mut buf,
            &header,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8)])))),
            b"8",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8), Some(13)])))),
            b"8,13",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8), None])))),
            b"8,.",
        )?;

        t(
            &mut buf,
            &header,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333)])))),
            b"0.333",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333), Some(0.667)])))),
            b"0.333,0.667",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333), None])))),
            b"0.333,.",
        )?;

        t(
            &mut buf,
            &header,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n')])))),
            b"n",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n'), Some('d')])))),
            b"n,d",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n'), None])))),
            b"n,.",
        )?;

        t(
            &mut buf,
            &header,
            &Value::Array(Array::String(Box::new(Values(&[Some(String::from(
                "noodles",
            ))])))),
            b"noodles",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::String(Box::new(Values(&[
                Some(String::from("noodles")),
                Some(String::from("vcf")),
            ])))),
            b"noodles,vcf",
        )?;
        t(
            &mut buf,
            &header,
            &Value::Array(Array::String(Box::new(Values(&[
                Some(String::from("noodles")),
                None,
            ])))),
            b"noodles,.",
        )?;

        Ok(())
    }
}
