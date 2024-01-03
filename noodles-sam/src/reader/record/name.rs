use std::{error, fmt};

use crate::alignment::record_buf::Name;

/// An error when a raw SAM record name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
        }
    }
}

pub(super) fn parse_name(src: &[u8], name: &mut Option<Name>) -> Result<(), ParseError> {
    if src.is_empty() {
        return Err(ParseError::Empty);
    }

    if let Some(name) = name {
        let buf = name.as_mut();
        buf.clear();
        buf.extend(src);
    } else {
        *name = Some(Name::from(src));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_name() -> Result<(), ParseError> {
        let mut name = None;

        parse_name(b"r0", &mut name)?;
        assert_eq!(name, Some(Name::from(b"r0")));

        assert_eq!(parse_name(b"", &mut name), Err(ParseError::Empty));

        Ok(())
    }
}
