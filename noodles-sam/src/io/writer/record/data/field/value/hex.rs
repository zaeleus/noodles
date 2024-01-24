use std::io::{self, Write};

pub(super) fn write_hex<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_hex() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_hex(&mut buf, b"CAFE")?;
        assert_eq!(buf, b"CAFE");

        buf.clear();
        write_hex(&mut buf, b"ndls")?;
        assert_eq!(buf, b"ndls");

        Ok(())
    }
}
