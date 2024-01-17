use std::io::{self, Write};

pub(super) fn write_tag<W>(writer: &mut W, tag: [u8; 2]) -> io::Result<()>
where
    W: Write,
{
    if is_valid(tag) {
        writer.write_all(&tag)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid tag"))
    }
}

// § 1.3 "The header section" (2023-05-24): "...`[A-Za-z][A-Za-z0-9]`..."
fn is_valid(tag: [u8; 2]) -> bool {
    tag[0].is_ascii_alphabetic() && tag[1].is_ascii_alphanumeric()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_tag() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_tag(&mut buf, [b'N', b'D'])?;
        assert_eq!(buf, b"ND");

        buf.clear();
        assert!(matches!(
            write_tag(&mut buf, [b'0', b'D']),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid([b'N', b'D']));
        assert!(is_valid([b'N', b'd']));
        assert!(is_valid([b'N', b'0']));
        assert!(is_valid([b'n', b'D']));
        assert!(is_valid([b'n', b'd']));
        assert!(is_valid([b'n', b'0']));

        assert!(!is_valid([b'0', b'D']));
        assert!(!is_valid([b'*', b'D']));
        assert!(!is_valid([b'D', b'*']));
    }
}
