mod field;

use std::io::{self, Write};

use noodles_gff::feature::record::Attributes;

use self::field::write_field;

pub(super) fn write_attributes<W>(writer: &mut W, attributes: &dyn Attributes) -> io::Result<()>
where
    W: Write,
{
    for (i, result) in attributes.iter().enumerate() {
        let (key, value) = result?;

        if i > 0 {
            write_separator(writer)?;
        }

        write_field(writer, &key, &value)?;
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b' ';
    writer.write_all(&[SEPARATOR])
}

#[cfg(test)]
mod tests {
    use bstr::BString;
    use noodles_gff::feature::record_buf::{Attributes as AttributesBuf, attributes::field::Value};

    use super::*;

    #[test]
    fn test_write_attributes() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, attributes: &AttributesBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_attributes(buf, attributes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let attributes = [(BString::from("gene_id"), Value::from("g0"))]
            .into_iter()
            .collect();
        t(&mut buf, &attributes, br#"gene_id "g0";"#)?;

        let attributes = [
            (BString::from("gene_id"), Value::from("g0")),
            (BString::from("gene_name"), Value::from("n0")),
        ]
        .into_iter()
        .collect();
        t(&mut buf, &attributes, br#"gene_id "g0"; gene_name "n0";"#)?;

        Ok(())
    }
}
