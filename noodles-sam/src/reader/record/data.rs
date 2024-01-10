pub(crate) mod field;

use std::{error, fmt};

use self::field::parse_field;
use crate::{alignment::record::data::field::Tag, record::Data};

/// An error when raw SAM record data fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::DuplicateTag(_) => None,
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag:?}"),
            Self::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(tag) = e.tag() {
                    write!(f, ": {tag:?}")?;
                }

                Ok(())
            }
        }
    }
}

pub(super) fn parse_data(mut src: &[u8], data: &mut Data) -> Result<(), ParseError> {
    while !src.is_empty() {
        let (tag, value) = parse_field(&mut src).map_err(ParseError::InvalidField)?;

        if let Some((t, _)) = data.insert(tag, value) {
            return Err(ParseError::DuplicateTag(t));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_data() -> Result<(), ParseError> {
        use crate::alignment::{record::data::field::tag, record_buf::data::field::Value};

        let mut data = Data::default();

        parse_data(b"", &mut data)?;
        assert!(data.is_empty());

        let nh = (tag::ALIGNMENT_HIT_COUNT, Value::from(1u8));
        let co = (tag::COMMENT, Value::from("ndls"));

        data.clear();
        parse_data(b"NH:i:1", &mut data)?;
        let expected = [nh.clone()].into_iter().collect();
        assert_eq!(data, expected);

        data.clear();
        parse_data(b"NH:i:1\tCO:Z:ndls", &mut data)?;
        let expected = [nh, co].into_iter().collect();
        assert_eq!(data, expected);

        data.clear();
        assert_eq!(
            parse_data(b"NH:i:1\tNH:i:1", &mut data),
            Err(ParseError::DuplicateTag(tag::ALIGNMENT_HIT_COUNT))
        );

        data.clear();
        assert!(matches!(
            parse_data(b"NH:i:ndls", &mut data),
            Err(ParseError::InvalidField(_))
        ));

        Ok(())
    }
}
