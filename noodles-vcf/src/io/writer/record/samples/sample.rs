use std::io::{self, Write};

use crate::{
    io::writer::record::MISSING,
    variant::record_buf::samples::{
        sample::{value::Array, Value},
        Sample,
    },
};

pub(super) fn write_sample<W>(writer: &mut W, sample: &Sample<'_>) -> io::Result<()>
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
    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => writer.write_all(s.as_bytes()),
        Value::Array(array) => write_array_value(writer, array),
    }
}

fn write_array_value<W>(writer: &mut W, array: &Array) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b",";

    match array {
        Array::Integer(values) => {
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
        }
        Array::Float(values) => {
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
        }
        Array::Character(values) => {
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
        }
        Array::String(values) => {
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
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::from(8), b"8")?;
        t(&mut buf, &Value::from(0.333), b"0.333")?;
        t(&mut buf, &Value::from('n'), b"n")?;
        t(&mut buf, &Value::from("noodles"), b"noodles")?;

        t(&mut buf, &Value::from(vec![Some(8)]), b"8")?;
        t(&mut buf, &Value::from(vec![Some(8), Some(13)]), b"8,13")?;
        t(&mut buf, &Value::from(vec![Some(8), None]), b"8,.")?;

        t(&mut buf, &Value::from(vec![Some(0.333)]), b"0.333")?;
        t(
            &mut buf,
            &Value::from(vec![Some(0.333), Some(0.667)]),
            b"0.333,0.667",
        )?;
        t(&mut buf, &Value::from(vec![Some(0.333), None]), b"0.333,.")?;

        t(&mut buf, &Value::from(vec![Some('n')]), b"n")?;
        t(&mut buf, &Value::from(vec![Some('n'), Some('d')]), b"n,d")?;
        t(&mut buf, &Value::from(vec![Some('n'), None]), b"n,.")?;

        t(
            &mut buf,
            &Value::from(vec![Some(String::from("noodles"))]),
            b"noodles",
        )?;
        t(
            &mut buf,
            &Value::from(vec![
                Some(String::from("noodles")),
                Some(String::from("vcf")),
            ]),
            b"noodles,vcf",
        )?;
        t(
            &mut buf,
            &Value::from(vec![Some(String::from("noodles")), None]),
            b"noodles,.",
        )?;

        Ok(())
    }
}
