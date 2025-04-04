use std::io::{self, Write};

use super::percent_encode;
use crate::feature::record::attributes::field::Value;

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    for (i, result) in value.iter().enumerate() {
        let v = result?;

        if i > 0 {
            write_separator(writer)?;
        }

        let s = percent_encode(v.as_ref());
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
    use bstr::BString;

    use super::*;
    use crate::feature::record_buf::attributes::field::Value as ValueBuf;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: &ValueBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, &value.into())?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &ValueBuf::from("noodles"), b"noodles")?;
        t(&mut buf, &ValueBuf::from("8,13"), b"8%2C13")?;
        t(
            &mut buf,
            &ValueBuf::from(vec![BString::from("8"), BString::from("13")]),
            b"8,13",
        )?;

        Ok(())
    }
}
