use std::io::{self, Write};

use bstr::BStr;

use super::percent_encode;

pub(super) fn write_type<W>(writer: &mut W, ty: &BStr) -> io::Result<()>
where
    W: Write,
{
    let s = percent_encode(ty);
    writer.write_all(s.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, ty: &BStr, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_type(buf, ty)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, BStr::new("exon"), b"exon")?;
        t(&mut buf, BStr::new("ex\ton"), b"ex%09on")?;
        t(&mut buf, BStr::new("ex\non"), b"ex%0Aon")?;
        t(&mut buf, BStr::new("100%"), b"100%25")?;
        t(&mut buf, BStr::new("SO:0000704"), b"SO:0000704")?;

        Ok(())
    }
}
