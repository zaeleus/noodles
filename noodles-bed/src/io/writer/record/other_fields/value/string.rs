use std::io::{self, Write};

use bstr::BStr;

pub(super) fn write_string<W>(writer: &mut W, s: &BStr) -> io::Result<()>
where
    W: Write,
{
    if is_valid(s) {
        writer.write_all(s)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid string",
        ))
    }
}

// § 1.3."Terminology and concepts" (2022-01-05): "All **fields** are 7-bit US ASCII printable
// characters."
fn is_valid(buf: &BStr) -> bool {
    fn is_valid_char(b: u8) -> bool {
        matches!(b, b' '..=b'~')
    }

    buf.iter().copied().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, b"ndls".as_bstr())?;
        assert_eq!(buf, b"ndls");

        buf.clear();
        assert!(matches!(
            write_string(&mut buf, b"\t".as_bstr()),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
