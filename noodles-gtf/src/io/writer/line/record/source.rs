use std::io::{self, Write};

use bstr::BStr;

pub(super) fn write_source<W>(writer: &mut W, source: &BStr) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(source)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_source() -> io::Result<()> {
        let mut buf = Vec::new();
        write_source(&mut buf, BStr::new("NDLS"))?;
        assert_eq!(buf, b"NDLS");
        Ok(())
    }
}
