use std::io::{self, Write};

use bstr::BStr;

use super::percent_encode;

pub(super) fn write_type<W>(writer: &mut W, ty: &BStr) -> io::Result<()>
where
    W: Write,
{
    let src = percent_encode(ty);
    writer.write_all(src.as_bytes())
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
        t(&mut buf, BStr::new("exon%"), b"exon%25")?;

        Ok(())
    }
}
