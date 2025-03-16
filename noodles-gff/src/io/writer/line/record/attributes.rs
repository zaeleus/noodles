mod field;

use std::io::{self, Write};

use self::field::write_field;
use super::write_missing;
use crate::feature::record::Attributes;

pub(super) fn write_attributes<W>(writer: &mut W, attributes: &dyn Attributes) -> io::Result<()>
where
    W: Write,
{
    if attributes.is_empty() {
        write_missing(writer)?;
    } else {
        for (i, result) in attributes.iter().enumerate() {
            let (tag, value) = result?;

            if i > 0 {
                write_separator(writer)?;
            }

            write_field(writer, tag.as_ref(), &value)?;
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
    use crate::feature::record_buf::{attributes::field::Value, Attributes as AttributesBuf};

    #[test]
    fn test_write_attributes() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, attributes: &AttributesBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_attributes(buf, attributes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &AttributesBuf::default(), b".")?;

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
