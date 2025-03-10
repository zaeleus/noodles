use std::io::{self, Write};

pub(super) fn write_comment<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    writer.write_all(s.as_bytes())
}

pub(super) fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const PREFIX: u8 = b'#';
    writer.write_all(&[PREFIX])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_comment() -> io::Result<()> {
        let mut buf = Vec::new();
        write_comment(&mut buf, "noodles")?;
        assert_eq!(buf, b"#noodles");
        Ok(())
    }
}
