use std::io::{self, Write};

use crate::variant::record::ReferenceBases;

pub(super) fn write_reference_bases<W, B>(writer: &mut W, reference_bases: B) -> io::Result<()>
where
    W: Write,
    B: ReferenceBases,
{
    for result in reference_bases.iter() {
        let base = result?;

        if is_valid(base) {
            writer.write_all(&[base])?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference base",
            ));
        }
    }

    Ok(())
}

// § 1.6.1.4 "Fixed fields: REF" (2023-08-23): "Each base must be one of A,C,G,T,N (case insensitive)."
fn is_valid(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_bases() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_reference_bases(&mut buf, "A")?;
        assert_eq!(buf, b"A");

        buf.clear();
        write_reference_bases(&mut buf, "AC")?;
        assert_eq!(buf, b"AC");

        buf.clear();
        assert!(matches!(
            write_reference_bases(&mut buf, "Z"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(b'A'));
        assert!(is_valid(b'C'));
        assert!(is_valid(b'G'));
        assert!(is_valid(b'T'));
        assert!(is_valid(b'N'));

        assert!(is_valid(b'a'));
        assert!(is_valid(b'c'));
        assert!(is_valid(b'g'));
        assert!(is_valid(b't'));
        assert!(is_valid(b'n'));

        assert!(!is_valid(b'Z'));
        assert!(!is_valid(b'z'));
    }
}
