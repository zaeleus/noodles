use std::io::{self, Write};

use bstr::BStr;

use super::percent_encode;

pub(super) fn write_source<W>(writer: &mut W, source: &BStr) -> io::Result<()>
where
    W: Write,
{
    let s = percent_encode(source);
    writer.write_all(s.as_bytes())
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
        t(&mut buf, BStr::new("NDLS\tv1"), b"NDLS%09v1")?;
        t(&mut buf, BStr::new("NDLS\nv1"), b"NDLS%0Av1")?;
        t(&mut buf, BStr::new("NDLS\rv1"), b"NDLS%0Dv1")?;
        t(&mut buf, BStr::new("NDLS\x7f"), b"NDLS%7F")?;
        t(&mut buf, BStr::new("50%"), b"50%25")?;

        // Neither spaces nor the characters reserved in column 9 are encoded.
        t(&mut buf, BStr::new("NDLS v1"), b"NDLS v1")?;
        t(&mut buf, BStr::new("a;b=c&d,e"), b"a;b=c&d,e")?;

        Ok(())
    }
}
