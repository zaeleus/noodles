mod field;

use std::{error, fmt};

use noodles_core as core;

use self::field::parse_field;
use crate::{
    record::{info::field::key, Info},
    Header,
};

/// An error when raw VCF record genotypes values fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A key is invalid.
    InvalidKey(key::ParseError),
    /// A value is missing.
    MissingValue,
    /// A value is invalid.
    InvalidValue(field::value::ParseError),
    /// A key is duplicated.
    DuplicateKey,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Empty => write!(f, "empty input"),
            ParseError::InvalidKey(_) => write!(f, "invalid key"),
            ParseError::MissingValue => write!(f, "missing value"),
            ParseError::InvalidValue(_) => write!(f, "invalid value"),
            ParseError::DuplicateKey => write!(f, "duplicate key"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_info(header: &Header, s: &str, info: &mut Info) -> Result<(), ParseError> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    for raw_field in s.split(DELIMITER) {
        let (key, value) = parse_field(header, raw_field)?;

        if info.insert(key, value).is_some() {
            return Err(ParseError::DuplicateKey);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info() -> Result<(), ParseError> {
        use crate::record::info::field::{key, Value};

        let header = Header::default();
        let mut info = Info::default();

        info.clear();
        parse_info(&header, "NS=2", &mut info)?;
        let expected = [(key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)))]
            .into_iter()
            .collect();
        assert_eq!(info, expected);

        info.clear();
        parse_info(&header, "NS=2;AA=T", &mut info)?;
        let expected = [
            (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2))),
            (
                key::ANCESTRAL_ALLELE,
                Some(Value::String(String::from("T"))),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(info, expected);

        info.clear();
        assert!(matches!(
            parse_info(&header, ".", &mut info),
            Err(ParseError::InvalidKey(_))
        ));

        info.clear();
        assert!(matches!(
            parse_info(&header, "NS=ndls", &mut info),
            Err(ParseError::InvalidValue(_))
        ));

        Ok(())
    }
}
