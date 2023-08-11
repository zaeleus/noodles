pub(crate) mod string;

use std::{borrow::Cow, error, fmt};

use self::string::{parse_escaped_string, parse_raw_string};

/// An error returned when a VCF header record map field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidString(string::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidString(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidString(_) => write!(f, "invalid string"),
        }
    }
}

pub fn parse_value<'a>(src: &mut &'a [u8]) -> Result<Cow<'a, str>, ParseError> {
    const QUOTATION_MARK: u8 = b'"';

    if let Some(buf) = src.strip_prefix(&[QUOTATION_MARK]) {
        *src = buf;

        parse_escaped_string(src).map_err(ParseError::InvalidString)
    } else {
        parse_raw_string(src)
            .map(Cow::from)
            .map_err(ParseError::InvalidString)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() {
        let mut src = &b"noodles-vcf,"[..];
        assert_eq!(parse_value(&mut src), Ok(Cow::from("noodles-vcf")));

        let mut src = &br#""noodles-vcf","#[..];
        assert_eq!(parse_value(&mut src), Ok(Cow::from("noodles-vcf")));

        let mut src = &b"noodles-vcf"[..];
        assert!(matches!(
            parse_value(&mut src),
            Err(ParseError::InvalidString(_))
        ));
    }
}
