use std::io::{self, Write};

use bstr::BStr;

use super::percent_encode;

pub(super) fn write_tag<W>(writer: &mut W, tag: &BStr) -> io::Result<()>
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
        fn t(buf: &mut Vec<u8>, key: &BStr, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_tag(buf, key)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, BStr::new("ID"), b"ID")?;
        t(&mut buf, BStr::new("%s"), b"%25s")?;

        Ok(())
    }
}
