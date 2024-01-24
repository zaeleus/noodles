use std::io::{self, Write};

pub(super) fn write_string<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if is_valid(buf) {
        writer.write_all(buf)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid string",
        ))
    }
}

// § 1.5 "The alignment section: optional fields" (2023-05-24): "`[ !-~]*`".
fn is_valid(buf: &[u8]) -> bool {
    fn is_valid_char(b: u8) -> bool {
        matches!(b, b' '..=b'~')
    }

    buf.iter().copied().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, b"ndls")?;
        assert_eq!(buf, b"ndls");

        buf.clear();
        assert!(matches!(
            write_string(&mut buf, b"nd\tls"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(b""));
        assert!(is_valid(b"ndls"));
        assert!(is_valid(b" "));

        assert!(!is_valid(b"\t"));
        assert!(!is_valid("🍜".as_bytes()));
    }
}
