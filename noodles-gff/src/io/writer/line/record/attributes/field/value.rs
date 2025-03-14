use std::io::{self, Write};

use super::percent_encode;
use crate::feature::record_buf::attributes::field::Value;

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    for (i, v) in value.iter().enumerate() {
        if i > 0 {
            write_separator(writer)?;
        }

        let s = percent_encode(v);
        writer.write_all(s.as_bytes())?;
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b',';
    writer.write_all(&[SEPARATOR])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, key: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, key)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::from("noodles"), b"noodles")?;
        t(&mut buf, &Value::from("8,13"), b"8%2C13")?;
        t(
            &mut buf,
            &Value::from(vec![String::from("8"), String::from("13")]),
            b"8,13",
        )?;

        Ok(())
    }
}
