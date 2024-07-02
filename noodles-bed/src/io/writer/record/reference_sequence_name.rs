use std::io::{self, Write};

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: &[u8],
) -> io::Result<()>
where
    W: Write,
{
    if is_valid(reference_sequence_name) {
        writer.write_all(reference_sequence_name)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid reference sequence name",
        ))
    }
}

// 1.5 "BED fields" (2022-01-05): `[[:alnum:]_]{1,255}`.
fn is_valid(buf: &[u8]) -> bool {
    fn is_valid_char(b: u8) -> bool {
        const UNDERSCORE: u8 = b'_';
        b.is_ascii_alphanumeric() || b == UNDERSCORE
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
        write_reference_sequence_name(&mut buf, b"sq0")?;
        assert_eq!(buf, b"sq0");

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, "🍜".as_bytes()),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
