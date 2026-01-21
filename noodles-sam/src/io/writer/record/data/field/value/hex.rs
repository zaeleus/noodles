use std::io::{self, Write};

pub(super) fn write_hex<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if is_valid(buf) {
        writer.write_all(buf)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid hex"))
    }
}

// § 1.5 "The alignment section: optional fields" (2024-11-06): "`([0-9A-F][0-9A-F])*`".
fn is_valid(buf: &[u8]) -> bool {
    fn is_even(n: usize) -> bool {
        n.is_multiple_of(2)
    }

    fn is_hexdigit(b: u8) -> bool {
        matches!(b, b'0'..=b'9' | b'A'..=b'F')
    }

    is_even(buf.len()) && buf.iter().copied().all(is_hexdigit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_hex() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_hex(&mut buf, b"CAFE")?;
        assert_eq!(buf, b"CAFE");

        buf.clear();
        assert!(matches!(
            write_hex(&mut buf, b"ndls"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(b""));
        assert!(is_valid(b"CAFE"));

        assert!(!is_valid(b" "));
        assert!(!is_valid(b"\t"));
        assert!(!is_valid(b"cafe"));
        assert!(!is_valid("🍜".as_bytes()));
    }
}
