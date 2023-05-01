use std::io::{self, Write};

use super::MISSING;
use crate::record::{
    genotypes::{
        sample::{value::Array, Value},
        Keys, Sample,
    },
    Genotypes,
};

pub(super) fn write_genotypes<W>(writer: &mut W, genotypes: &Genotypes) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b"\t";

    write_keys(writer, genotypes.keys())?;

    for sample in genotypes.values() {
        writer.write_all(DELIMITER)?;
        write_sample(writer, sample)?;
    }

    Ok(())
}

fn write_keys<W>(writer: &mut W, keys: &Keys) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b":";

    for (i, key) in keys.iter().enumerate() {
        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        writer.write_all(key.as_ref().as_bytes())?;
    }

    Ok(())
}

fn write_sample<W>(writer: &mut W, sample: Sample<'_>) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b":";

    for (i, value) in sample.values().iter().enumerate() {
        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        match value {
            Some(v) => write_value(writer, v)?,
            None => writer.write_all(MISSING)?,
        }
    }

    Ok(())
}

fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b",";

    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => writer.write_all(s.as_bytes()),
        Value::Array(Array::Integer(values)) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = v {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::Array(Array::Float(values)) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = v {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::Array(Array::Character(values)) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(c) = v {
                    write!(writer, "{c}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::Array(Array::String(values)) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(s) = v {
                    writer.write_all(s.as_bytes())?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::genotypes::keys::key;

        fn t(buf: &mut Vec<u8>, genotypes: &Genotypes, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_genotypes(buf, genotypes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let genotypes = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE])?,
            vec![vec![Some(Value::String(String::from("0|0")))]],
        );
        t(&mut buf, &genotypes, b"GT\t0|0")?;

        let genotypes = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![
                vec![
                    Some(Value::String(String::from("0|0"))),
                    Some(Value::Integer(13)),
                ],
                vec![
                    Some(Value::String(String::from("0/1"))),
                    Some(Value::Integer(8)),
                ],
            ],
        );
        t(&mut buf, &genotypes, b"GT:GQ\t0|0:13\t0/1:8")?;

        Ok(())
    }

    #[test]
    fn test_write_value() -> io::Result<()> {
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
        t(
            &mut buf,
            &Value::String(String::from("noodles")),
            b"noodles",
        )?;

        t(&mut buf, &Value::Array(Array::Integer(vec![Some(8)])), b"8")?;
        t(
            &mut buf,
            &Value::Array(Array::Integer(vec![Some(8), Some(13)])),
            b"8,13",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Integer(vec![Some(8), None])),
            b"8,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::Float(vec![Some(0.333)])),
            b"0.333",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Float(vec![Some(0.333), Some(0.667)])),
            b"0.333,0.667",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Float(vec![Some(0.333), None])),
            b"0.333,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::Character(vec![Some('n')])),
            b"n",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Character(vec![Some('n'), Some('d')])),
            b"n,d",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::Character(vec![Some('n'), None])),
            b"n,.",
        )?;

        t(
            &mut buf,
            &Value::Array(Array::String(vec![Some(String::from("noodles"))])),
            b"noodles",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::String(vec![
                Some(String::from("noodles")),
                Some(String::from("vcf")),
            ])),
            b"noodles,vcf",
        )?;
        t(
            &mut buf,
            &Value::Array(Array::String(vec![Some(String::from("noodles")), None])),
            b"noodles,.",
        )?;

        Ok(())
    }
}
