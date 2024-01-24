use std::io::{self, Write};

pub(super) fn write_string<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, b"ndls")?;
        assert_eq!(buf, b"ndls");

        buf.clear();
        write_string(&mut buf, b"nd\tls")?;
        assert_eq!(buf, b"nd\tls");

        Ok(())
    }
}
