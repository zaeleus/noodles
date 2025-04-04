use std::io::{self, Write};

use bstr::BStr;

pub(super) fn write_comment<W>(writer: &mut W, s: &BStr) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    writer.write_all(s)
}

pub(super) fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const PREFIX: u8 = b'#';
    writer.write_all(&[PREFIX])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_comment() -> io::Result<()> {
        let mut buf = Vec::new();
        write_comment(&mut buf, BStr::new("noodles"))?;
        assert_eq!(buf, b"#noodles");
        Ok(())
    }
}
