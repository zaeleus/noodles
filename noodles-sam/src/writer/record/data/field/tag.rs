use std::io::{self, Write};

pub fn write_tag<W>(writer: &mut W, tag: [u8; 2]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&tag[..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() -> io::Result<()> {
        let mut buf = Vec::new();
        write_tag(&mut buf, [b'N', b'H'])?;
        assert_eq!(buf, b"NH");
        Ok(())
    }
}
