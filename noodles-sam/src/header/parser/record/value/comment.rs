use std::{error, fmt, str};

/// An error returned when a SAM header record comment value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_comment(src: &mut &[u8]) -> Result<String, ParseError> {
    let s = str::from_utf8(src)
        .map(String::from)
        .map_err(ParseError::Invalid)?;

    *src = &[];

    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_comment() {
        let mut src = &b"noodles"[..];
        assert_eq!(parse_comment(&mut src), Ok(String::from("noodles")));

        let mut src = &[0x00, 0x9f, 0x8d, 0x9c][..];
        assert!(matches!(
            parse_comment(&mut src),
            Err(ParseError::Invalid(_))
        ));
    }
}
