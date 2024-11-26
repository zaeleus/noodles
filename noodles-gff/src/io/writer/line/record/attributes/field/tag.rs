use std::io::{self, Write};

use super::percent_encode;

pub(super) fn write_tag<W>(writer: &mut W, tag: &str) -> io::Result<()>
where
    W: Write,
{
    let s = percent_encode(tag);
    writer.write_all(s.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_tag() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, key: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_tag(buf, key)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "ID", b"ID")?;
        t(&mut buf, "%s", b"%25s")?;

        Ok(())
    }
}
