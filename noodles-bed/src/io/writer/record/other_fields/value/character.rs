use std::io::{self, Write};

pub(super) fn write_character<W>(writer: &mut W, b: u8) -> io::Result<()>
where
    W: Write,
{
    if is_valid(b) {
        writer.write_all(&[b])
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid character",
        ))
    }
}

// § 1.3."Terminology and concepts" (2022-01-05): "All **fields** are 7-bit US ASCII printable
// characters."
fn is_valid(b: u8) -> bool {
    matches!(b, b' '..=b'~')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_character() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_character(&mut buf, b'n')?;
        assert_eq!(buf, b"n");

        buf.clear();
        assert!(matches!(
            write_character(&mut buf, b'\t'),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
