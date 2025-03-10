use std::io::{self, Write};

use super::write_separator;
use crate::record_buf::attributes::Entry;

pub(super) fn write_field<W>(writer: &mut W, field: &Entry) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(field.key().as_bytes())?;
    write_separator(writer)?;
    write!(writer, r#""{}""#, field.value())?;
    write_terminator(writer)?;
    Ok(())
}

fn write_terminator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const TERMINATOR: u8 = b';';
    writer.write_all(&[TERMINATOR])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let field = Entry::new("gene_id", "g0");
        write_field(&mut buf, &field)?;
        assert_eq!(buf, br#"gene_id "g0";"#);
        Ok(())
    }
}
