use std::io::{self, Write};

use crate::{
    io::writer::record::MISSING,
    variant::record::samples::{
        series::{
            value::{genotype::Phasing, Array, Genotype},
            Value,
        },
        Sample,
    },
    Header,
};

pub(super) fn write_sample<W, S>(writer: &mut W, header: &Header, sample: S) -> io::Result<()>
where
    W: Write,
    S: Sample,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in sample.iter(header).enumerate() {
        let (_, value) = result?;

        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        match value {
            Some(v) => write_value(writer, &v)?,
            None => writer.write_all(MISSING)?,
        }
    }

    Ok(())
}

fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => writer.write_all(s.as_bytes()),
        Value::Genotype(genotype) => write_genotype(writer, genotype.as_ref()),
        Value::Array(array) => write_array_value(writer, array),
    }
}

fn write_genotype<W>(writer: &mut W, genotype: &dyn Genotype) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = b'.';

    for (i, result) in genotype.iter().enumerate() {
        let (position, phasing) = result?;

        if i > 0 {
            write_genotype_phasing(writer, phasing)?;
        }

        if let Some(n) = position {
            write!(writer, "{n}")?;
        } else {
            writer.write_all(&[MISSING])?;
        }
    }

    Ok(())
}

fn write_genotype_phasing<W>(writer: &mut W, phasing: Phasing) -> io::Result<()>
where
    W: Write,
{
    const PHASED: u8 = b'/';
    const UNPHASED: u8 = b'|';

    match phasing {
        Phasing::Phased => writer.write_all(&[PHASED]),
        Phasing::Unphased => writer.write_all(&[UNPHASED]),
    }
}

fn write_array_value<W>(writer: &mut W, array: &Array) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b",";

    match array {
        Array::Integer(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = result? {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::Float(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = result? {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::Character(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(c) = result? {
                    write!(writer, "{c}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::String(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(s) = result? {
                    writer.write_all(s.as_bytes())?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Integer(8), b"8")?;
        t(&mut buf, &Value::Float(0.333), b"0.333")?;
        t(&mut buf, &Value::Character('n'), b"n")?;
        t(&mut buf, &Value::String("noodles"), b"noodles")?;

        t(
            &mut buf,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8)])))),
            b"8",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8), Some(13)])))),
            b"8,13",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Integer(Box::new(Values(&[Some(8), None])))),
            b"8,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333)])))),
            b"0.333",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333), Some(0.667)])))),
            b"0.333,0.667",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Float(Box::new(Values(&[Some(0.333), None])))),
            b"0.333,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n')])))),
            b"n",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n'), Some('d')])))),
            b"n,d",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Character(Box::new(Values(&[Some('n'), None])))),
            b"n,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::String(Box::new(Values(&[Some(String::from(
                "noodles",
            ))])))),
            b"noodles",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::String(Box::new(Values(&[
                Some(String::from("noodles")),
                Some(String::from("vcf")),
            ])))),
            b"noodles,vcf",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::String(Box::new(Values(&[
                Some(String::from("noodles")),
                None,
            ])))),
            b"noodles,.",
        )?;

        Ok(())
    }
}
