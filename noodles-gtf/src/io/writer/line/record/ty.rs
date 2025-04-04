use std::io::{self, Write};

use bstr::BStr;

pub(super) fn write_type<W>(writer: &mut W, ty: &BStr) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(ty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, BStr::new("exon"))?;
        assert_eq!(buf, b"exon");
        Ok(())
    }
}
