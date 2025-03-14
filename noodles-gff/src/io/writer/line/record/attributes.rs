mod field;

use std::io::{self, Write};

use self::field::write_field;
use super::write_missing;
use crate::feature::record_buf::Attributes;

pub(super) fn write_attributes<W>(writer: &mut W, attributes: &Attributes) -> io::Result<()>
where
    W: Write,
{
    if attributes.is_empty() {
        write_missing(writer)?;
    } else {
        for (i, (tag, value)) in attributes.as_ref().iter().enumerate() {
            if i > 0 {
                write_separator(writer)?;
            }

            write_field(writer, tag, value)?;
        }
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b';';
    writer.write_all(&[SEPARATOR])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::record_buf::attributes::field::Value;

    #[test]
    fn test_write_attributes() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, attributes: &Attributes, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_attributes(buf, attributes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Attributes::default(), b".")?;

        t(
            &mut buf,
            &[(String::from("ID"), Value::from("0"))]
                .into_iter()
                .collect(),
            b"ID=0",
        )?;

        t(
            &mut buf,
            &[
                (String::from("ID"), Value::from("0")),
                (String::from("Name"), Value::from("ndls")),
            ]
            .into_iter()
            .collect(),
            b"ID=0;Name=ndls",
        )?;

        Ok(())
    }
}
