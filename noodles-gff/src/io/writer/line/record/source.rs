use std::io::{self, Write};

use bstr::BStr;

use super::percent_encode;

pub(super) fn write_source<W>(writer: &mut W, source: &BStr) -> io::Result<()>
where
    W: Write,
{
    let src = percent_encode(source);
    writer.write_all(src.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_source() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, source: &BStr, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_source(buf, source)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, BStr::new("NDLS"), b"NDLS")?;
        t(&mut buf, BStr::new("NDLS%"), b"NDLS%25")?;

        Ok(())
    }
}
