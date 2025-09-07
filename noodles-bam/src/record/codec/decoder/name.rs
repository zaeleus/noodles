use std::{
    error, fmt,
    num::{self, NonZero},
};

use bstr::BString;

const NUL: u8 = 0x00;

/// An error when a raw BAM record name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
    /// The NUL terminator is missing.
    MissingNulTerminator { actual: u8 },
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidLength(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::MissingNulTerminator { actual } => write!(
                f,
                "missing NUL terminator: expected {NUL:#04x}, got {actual:#04x}"
            ),
        }
    }
}

pub(super) fn read_length(src: &mut &[u8]) -> Result<NonZero<usize>, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    let len = usize::from(*n);
    *src = rest;
    NonZero::try_from(len).map_err(DecodeError::InvalidLength)
}

pub(super) fn read_name(
    src: &mut &[u8],
    name: &mut Option<BString>,
    len: NonZero<usize>,
) -> Result<(), DecodeError> {
    const MISSING: [u8; 2] = [b'*', NUL];

    let (buf, rest) = src
        .split_at_checked(len.get())
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    *name = if buf == MISSING {
        None
    } else {
        let mut dst = name.take().unwrap_or_default();

        // SAFETY: `buf` is non-empty.
        let (terminator, buf) = buf.split_last().unwrap();

        dst.resize(buf.len(), 0);
        dst.copy_from_slice(buf);

        if *terminator != NUL {
            return Err(DecodeError::MissingNulTerminator {
                actual: *terminator,
            });
        }

        Some(dst)
    };

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_length() {
        let mut src = &8i32.to_le_bytes()[..];
        assert_eq!(
            read_length(&mut src),
            Ok(const { NonZero::new(8).unwrap() })
        );

        let mut src = &[][..];
        assert_eq!(read_length(&mut src), Err(DecodeError::UnexpectedEof));

        let mut src = &0i32.to_le_bytes()[..];
        assert!(matches!(
            read_length(&mut src),
            Err(DecodeError::InvalidLength(_))
        ));
    }

    #[test]
    fn test_read_name() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], expected: Option<BString>) -> Result<(), DecodeError> {
            let mut actual = None;
            let l_read_name = NonZero::try_from(src.len()).unwrap();
            read_name(&mut src, &mut actual, l_read_name)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'*', 0x00], None)?;
        t(&[b'r', b'1', 0x00], Some(BString::from(b"r1")))?;

        let src = [0xf0, 0x9f, 0x8d, 0x9c, 0x00]; // "🍜\x00"
        t(&src, Some(BString::from(&src[0..4])))?;

        let data = [b'*'];
        let mut src = &data[..];
        let l_read_name = NonZero::try_from(data.len()).unwrap();
        assert_eq!(
            read_name(&mut src, &mut None, l_read_name),
            Err(DecodeError::MissingNulTerminator { actual: b'*' })
        );

        Ok(())
    }
}
