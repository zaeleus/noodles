use std::{error, fmt};

/// An error when raw SAM record flags fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(lexical_core::Error),
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

pub(crate) fn parse_template_length(src: &[u8]) -> Result<i32, ParseError> {
    lexical_core::parse(src).map_err(ParseError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_template_length() {
        assert_eq!(parse_template_length(b"0"), Ok(0));
        assert_eq!(parse_template_length(b"8"), Ok(8));

        assert!(matches!(
            parse_template_length(b""),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            parse_template_length(b"n"),
            Err(ParseError::Invalid(_))
        ));
    }
}
