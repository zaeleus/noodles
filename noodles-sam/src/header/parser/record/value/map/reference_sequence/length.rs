use std::{error, fmt, num::NonZeroUsize};

/// An error returned when a SAM header reference sequence record length value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(lexical_core::Error),
    /// The value is 0.
    Zero,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            Self::Zero => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
            Self::Zero => write!(f, "invalid value: 0"),
        }
    }
}

pub(super) fn parse_length(src: &mut &[u8]) -> Result<NonZeroUsize, ParseError> {
    let (n, i) = lexical_core::parse_partial::<usize>(src).map_err(ParseError::Invalid)?;
    *src = &src[i..];
    NonZeroUsize::new(n).ok_or(ParseError::Zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_length() {
        const LN_8: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let mut src = &b"8"[..];
        assert_eq!(parse_length(&mut src), Ok(LN_8));

        let mut src = &b"538522340430300790495419781092981030533"[..];
        assert!(matches!(
            parse_length(&mut src),
            Err(ParseError::Invalid(_))
        ));

        let mut src = &b"0"[..];
        assert_eq!(parse_length(&mut src), Err(ParseError::Zero));
    }
}
