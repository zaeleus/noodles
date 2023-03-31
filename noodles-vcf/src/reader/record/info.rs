use std::{error, fmt};

use noodles_core as core;

use crate::{
    record::{info::field, Info},
    Header,
};

/// An error when raw VCF record genotypes values fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is invalid.
    InvalidField,
    /// A key is duplicated.
    DuplicateKey,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Empty => write!(f, "empty input"),
            ParseError::InvalidField => write!(f, "invalid field"),
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
        let (key, value) =
            field::parse(raw_field, header.infos()).map_err(|_| ParseError::InvalidField)?;

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
        use crate::{header::info::key, record::info::field::Value};

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
        assert_eq!(
            parse_info(&header, ".", &mut info),
            Err(ParseError::InvalidField)
        );

        info.clear();
        assert_eq!(
            parse_info(&header, "NS=ndls", &mut info),
            Err(ParseError::InvalidField)
        );

        Ok(())
    }
}
