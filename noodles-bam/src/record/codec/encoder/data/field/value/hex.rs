use std::io;

use bstr::BStr;

use crate::record::codec::encoder::num::write_u8;

// § 4.2.4 " Auxiliary data encoding" (2024-11-06): "...hex-formatted byte arrays are
// represented as `NUL`-terminated text strings...".
pub(super) fn write_hex(dst: &mut Vec<u8>, s: &BStr) -> io::Result<()> {
    const NUL: u8 = 0x00;

    if is_valid(s) {
        dst.extend(s.iter());
        write_u8(dst, NUL);
        Ok(())
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid hex"))
    }
}

fn is_valid(s: &BStr) -> bool {
    fn is_even(n: usize) -> bool {
        n % 2 == 0
    }

    // § 1.5 "The alignment section: optional fields" (2024-11-06): "`([0-9A-F][0-9A-F])*`".
    fn is_hexdigit(b: u8) -> bool {
        matches!(b, b'0'..=b'9' | b'A'..=b'F')
    }

    is_even(s.len()) && s.iter().copied().all(is_hexdigit)
}

#[cfg(test)]
mod tests {
    use bstr::BStr;

    use super::*;

    #[test]
    fn test_write_hex() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_hex(&mut buf, BStr::new(b"CAFE"))?;
        assert_eq!(buf, [b'C', b'A', b'F', b'E', 0x00]);

        buf.clear();
        assert!(matches!(
            write_hex(&mut buf, BStr::new(b"ndls")),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(BStr::new(b"")));
        assert!(is_valid(BStr::new(b"CAFE")));

        assert!(!is_valid(BStr::new(b" ")));
        assert!(!is_valid(BStr::new(b"\t")));
        assert!(!is_valid(BStr::new(b"cafe")));
        assert!(!is_valid(BStr::new("🍜")));
    }
}
