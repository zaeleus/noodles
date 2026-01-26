use std::{error, fmt, mem};

use bstr::{BStr, BString};

use super::num::write_u8;

const MIN_LENGTH: usize = 1;
const MAX_LENGTH: usize = 254;
const MISSING: &[u8] = b"*";

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EncodeError {
    /// The name length is invalid.
    InvalidLength(usize),
    /// The name is invalid.
    Invalid(BString),
}

impl error::Error for EncodeError {}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(n) => write!(
                f,
                "invalid length: expected {MIN_LENGTH} <= n <= {MAX_LENGTH}, got {n}"
            ),
            Self::Invalid(buf) => write!(f, "invalid name: {buf}"),
        }
    }
}

pub(super) fn write_length(dst: &mut Vec<u8>, name: Option<&BStr>) -> Result<(), EncodeError> {
    let mut len = name.map(|s| s.len()).unwrap_or(MISSING.len());

    // + NUL terminator
    len += mem::size_of::<u8>();

    let n = u8::try_from(len).map_err(|_| EncodeError::InvalidLength(len))?;
    write_u8(dst, n);

    Ok(())
}

pub(super) fn write_name(dst: &mut Vec<u8>, name: Option<&BStr>) -> Result<(), EncodeError> {
    const NUL: u8 = 0x00;

    if let Some(name) = name {
        if !is_valid(name) {
            return Err(EncodeError::Invalid(name.into()));
        }

        dst.extend_from_slice(name.as_ref());
    } else {
        dst.extend(MISSING);
    }

    write_u8(dst, NUL);

    Ok(())
}

fn is_valid(buf: &[u8]) -> bool {
    (MIN_LENGTH..=MAX_LENGTH).contains(&buf.len())
        && buf != MISSING
        && buf.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_length() -> Result<(), EncodeError> {
        let mut buf = Vec::new();

        buf.clear();
        write_length(&mut buf, None)?;
        assert_eq!(buf, [0x02]);

        buf.clear();
        write_length(&mut buf, Some(b"r0".as_bstr()))?;
        assert_eq!(buf, [0x03]);

        buf.clear();
        let name = vec![b'n'; 255];
        assert!(matches!(
            write_length(&mut buf, Some(name.as_bstr())),
            Err(EncodeError::InvalidLength(256))
        ));

        Ok(())
    }

    #[test]
    fn test_write_name() -> Result<(), EncodeError> {
        fn t(buf: &mut Vec<u8>, name: Option<&[u8]>, expected: &[u8]) -> Result<(), EncodeError> {
            buf.clear();
            write_name(buf, name.map(|buf| buf.as_bstr()))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[b'*', 0x00])?;
        t(&mut buf, Some(b"r0"), &[b'r', b'0', 0x00])?;

        buf.clear();
        assert!(matches!(
            write_name(&mut buf, Some(MISSING.as_bstr())),
            Err(EncodeError::Invalid(_))
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(b"r0"));

        assert!(!is_valid(b""));
        assert!(!is_valid(b"*"));
        assert!(!is_valid(b"r 0"));
        assert!(!is_valid(b"@r0"));

        let s = vec![b'n'; MAX_LENGTH + 1];
        assert!(!is_valid(&s));
    }
}
