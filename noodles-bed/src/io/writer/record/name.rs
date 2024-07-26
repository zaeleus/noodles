use std::io::{self, Write};

pub(super) fn write_name<W>(writer: &mut W, name: Option<&[u8]>) -> io::Result<()>
where
    W: Write,
{
    const MISSING: &[u8] = b".";

    if let Some(buf) = name {
        if is_valid(buf) {
            writer.write_all(buf)
        } else {
            Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid name"))
        }
    } else {
        writer.write_all(MISSING)
    }
}

// 1.5 "BED fields" (2022-01-05): `[\x20-\x7e]{1,255}`.
fn is_valid(buf: &[u8]) -> bool {
    fn is_valid_char(b: u8) -> bool {
        matches!(b, b' '..=b'~')
    }

    matches!(buf.len(), 1..=255) && buf.iter().copied().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_name(&mut buf, None)?;
        assert_eq!(buf, b".");

        buf.clear();
        write_name(&mut buf, Some(b"r0"))?;
        assert_eq!(buf, b"r0");

        buf.clear();
        assert!(matches!(
            write_name(&mut buf, Some("🍜".as_bytes())),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
