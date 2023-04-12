pub(crate) mod field;

use std::{error, fmt};

use self::field::parse_field;
use crate::record::Data;

/// An error when raw SAM record data fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(tag) = e.tag() {
                    write!(f, ": {tag}")?;
                }
            }
        }

        Ok(())
    }
}

pub(crate) fn parse_data(mut src: &[u8], data: &mut Data) -> Result<(), ParseError> {
    while !src.is_empty() {
        let (tag, value) = parse_field(&mut src).map_err(ParseError::InvalidField)?;
        data.insert(tag, value);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_data() -> Result<(), ParseError> {
        use crate::record::data::field::{Tag, Value};

        let mut data = Data::default();

        parse_data(b"", &mut data)?;
        assert!(data.is_empty());

        let nh = (Tag::AlignmentHitCount, Value::from(1u8));
        let co = (Tag::Comment, Value::String(String::from("ndls")));

        data.clear();
        parse_data(b"NH:i:1", &mut data)?;
        let expected = [nh.clone()].into_iter().collect();
        assert_eq!(data, expected);

        data.clear();
        parse_data(b"NH:i:1\tCO:Z:ndls", &mut data)?;
        let expected = [nh, co].into_iter().collect();
        assert_eq!(data, expected);

        Ok(())
    }
}
