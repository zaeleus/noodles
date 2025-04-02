use std::io;

use bstr::BStr;

use crate::record::codec::encoder::num::write_u8;

// § 4.2.4 " Auxiliary data encoding" (2024-11-06): "String fields ... are represented as
// `NUL`-terminated text strings..."
pub(super) fn write_string(dst: &mut Vec<u8>, s: &BStr) -> io::Result<()> {
    const NUL: u8 = 0x00;

    if is_valid(s) {
        dst.extend(s.iter());
        write_u8(dst, NUL);
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid string",
        ))
    }
}

// § 1.5 "The alignment section: optional fields" (2024-11-06): "`[ !-~]*`".
fn is_valid(s: &BStr) -> bool {
    fn is_valid_char(b: u8) -> bool {
        matches!(b, b' '..=b'~')
    }

    s.iter().copied().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use bstr::BStr;

    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, BStr::new(b"ndls"))?;
        assert_eq!(buf, [b'n', b'd', b'l', b's', 0x00]);

        buf.clear();
        assert!(matches!(
            write_string(&mut buf, BStr::new("🍜")),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(BStr::new(b"")));
        assert!(is_valid(BStr::new(b"ndls")));
        assert!(is_valid(BStr::new(b" ")));

        assert!(!is_valid(BStr::new(b"\t")));
        assert!(!is_valid(BStr::new("🍜")));
    }
}
