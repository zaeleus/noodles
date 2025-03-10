use std::io::{self, Write};

use super::write_separator;

pub(super) fn write_field<W>(writer: &mut W, key: &str, value: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(key.as_bytes())?;
    write_separator(writer)?;
    write!(writer, r#""{}""#, value)?;
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
        write_field(&mut buf, "gene_id", "g0")?;
        assert_eq!(buf, br#"gene_id "g0";"#);
        Ok(())
    }
}
