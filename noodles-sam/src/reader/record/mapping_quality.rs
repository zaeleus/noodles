use std::{error, fmt};

use crate::alignment::record_buf::MappingQuality;

/// An error when a raw SAM record mapping quality fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
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

pub(super) fn parse_mapping_quality(src: &[u8]) -> Result<Option<MappingQuality>, ParseError> {
    lexical_core::parse(src)
        .map_err(ParseError::Invalid)
        .map(MappingQuality::new)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapping_quality() {
        assert_eq!(parse_mapping_quality(b"0"), Ok(MappingQuality::new(0)));
        assert_eq!(parse_mapping_quality(b"8"), Ok(MappingQuality::new(8)));
        assert_eq!(parse_mapping_quality(b"255"), Ok(None));

        assert!(matches!(
            parse_mapping_quality(b""),
            Err(ParseError::Invalid(_))
        ));
        assert!(matches!(
            parse_mapping_quality(b"256"),
            Err(ParseError::Invalid(_))
        ));
    }
}
