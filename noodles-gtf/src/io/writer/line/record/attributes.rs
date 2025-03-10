mod field;

use std::io::{self, Write};

use self::field::write_field;
use crate::record::Attributes;

pub(super) fn write_attributes<W>(writer: &mut W, attributes: &Attributes) -> io::Result<()>
where
    W: Write,
{
    for (i, field) in attributes.as_ref().iter().enumerate() {
        if i > 0 {
            write_separator(writer)?;
        }

        write_field(writer, field)?;
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
    use super::*;
    use crate::record::attributes::Entry;

    #[test]
    fn test_write_attributes() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, attributes: &Attributes, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_attributes(buf, attributes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let attributes = Attributes::from(vec![Entry::new("gene_id", "g0")]);
        t(&mut buf, &attributes, br#"gene_id "g0";"#)?;

        let attributes = Attributes::from(vec![
            Entry::new("gene_id", "g0"),
            Entry::new("gene_name", "n0"),
        ]);
        t(&mut buf, &attributes, br#"gene_id "g0"; gene_name "n0";"#)?;

        Ok(())
    }
}
