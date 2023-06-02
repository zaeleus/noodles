use std::{
    error, fmt, mem,
    num::{self, NonZeroUsize},
};

use bytes::Buf;
use noodles_sam::record::{read_name, ReadName};

const NUL: u8 = 0x00;

/// An error when a raw BAM record read name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
    /// The NUL terminator is missing.
    MissingNulTerminator { actual: u8 },
    /// An input is invalid.
    Invalid(read_name::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidLength(e) => Some(e),
            Self::Invalid(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::MissingNulTerminator { actual } => write!(
                f,
                "missing NUL terminator: expected {NUL:#04x}, got {actual:#04x}"
            ),
            Self::Invalid(_) => todo!(),
        }
    }
}

pub(crate) fn get_length<B>(src: &mut B) -> Result<NonZeroUsize, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(ParseError::UnexpectedEof);
    }

    NonZeroUsize::try_from(usize::from(src.get_u8())).map_err(ParseError::InvalidLength)
}

pub fn get_read_name<B>(
    src: &mut B,
    read_name: &mut Option<ReadName>,
    l_read_name: NonZeroUsize,
) -> Result<(), ParseError>
where
    B: Buf,
{
    const MISSING: [u8; 2] = [b'*', NUL];

    let len = usize::from(l_read_name);

    if src.remaining() < len {
        return Err(ParseError::UnexpectedEof);
    }

    *read_name = if src.take(len).chunk() == MISSING {
        src.advance(MISSING.len());
        None
    } else {
        let mut dst = read_name.take().map(Vec::from).unwrap_or_default();

        // SAFETY: len is guaranteed to be > 0.
        dst.resize(len - 1, 0);
        src.copy_to_slice(&mut dst);

        let terminator = src.get_u8();

        if terminator != NUL {
            return Err(ParseError::MissingNulTerminator { actual: terminator });
        }

        ReadName::try_from(dst)
            .map(Some)
            .map_err(ParseError::Invalid)?
    };

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_length() -> Result<(), num::TryFromIntError> {
        let mut src = &8i32.to_le_bytes()[..];
        assert_eq!(get_length(&mut src), Ok(NonZeroUsize::try_from(8)?));

        let mut src = &[][..];
        assert_eq!(get_length(&mut src), Err(ParseError::UnexpectedEof));

        let mut src = &0i32.to_le_bytes()[..];
        assert!(matches!(
            get_length(&mut src),
            Err(ParseError::InvalidLength(_))
        ));

        Ok(())
    }

    #[test]
    fn test_get_read_name() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], expected: Option<ReadName>) -> Result<(), ParseError> {
            let mut actual = None;
            let l_read_name = NonZeroUsize::try_from(src.len()).unwrap();
            get_read_name(&mut src, &mut actual, l_read_name)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'*', 0x00], None)?;
        t(&[b'r', 0x00], "r".parse().map(Some)?)?;
        t(&[b'r', b'1', 0x00], "r1".parse().map(Some)?)?;

        let data = [b'*'];
        let mut src = &data[..];
        let l_read_name = NonZeroUsize::try_from(data.len()).unwrap();
        assert_eq!(
            get_read_name(&mut src, &mut None, l_read_name),
            Err(ParseError::MissingNulTerminator { actual: b'*' })
        );

        let data = [0xf0, 0x9f, 0x8d, 0x9c, 0x00]; // "🍜\x00"
        let mut src = &data[..];
        let l_read_name = NonZeroUsize::try_from(data.len()).unwrap();
        assert!(matches!(
            get_read_name(&mut src, &mut None, l_read_name),
            Err(ParseError::Invalid(_))
        ));

        Ok(())
    }
}
