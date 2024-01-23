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

// § 1.5 "The alignment section: optional fields" (2023-05-24): "`[!-~]`".
fn is_valid(b: u8) -> bool {
    b.is_ascii_graphic()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_character() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_character(&mut buf, b'n')?;
        assert_eq!(buf, [b'n']);

        buf.clear();
        assert!(matches!(
            write_character(&mut buf, b'\n'),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid_character() {
        assert!(is_valid(b'n'));
        assert!(!is_valid(b'\t'));
    }
}
