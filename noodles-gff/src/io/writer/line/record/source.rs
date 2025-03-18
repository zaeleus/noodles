use std::io::{self, Write};

pub(super) fn write_source<W>(writer: &mut W, source: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(source.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_source() -> io::Result<()> {
        let mut buf = Vec::new();
        write_source(&mut buf, "NDLS")?;
        assert_eq!(buf, b"NDLS");
        Ok(())
    }
}
