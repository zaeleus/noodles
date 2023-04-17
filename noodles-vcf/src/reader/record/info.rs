mod value;

use std::{error, fmt};

use noodles_core as core;

use self::value::parse_value;
use crate::{
    header::{record::value::map::info::Type, Number},
    record::{
        info::field::{key, Key, Value},
        Info,
    },
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
    InvalidValue(value::ParseError),
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

fn parse_field(header: &Header, s: &str) -> Result<(Key, Option<Value>), ParseError> {
    use super::MISSING;

    const MAX_COMPONENTS: usize = 2;
    const SEPARATOR: char = '=';

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let raw_key = components.next().unwrap_or_default();
    let key = raw_key.parse().map_err(ParseError::InvalidKey)?;

    let (number, ty) = header
        .infos()
        .get(&key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| {
            crate::header::info::key::definition(header.file_format(), &key).map(|(n, t, _)| (n, t))
        })
        .unwrap_or((Number::Count(1), Type::String));

    let raw_value = components.next();

    let value = if matches!(ty, Type::Flag) {
        match raw_value.unwrap_or_default() {
            MISSING => None,
            t => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
        }
    } else if matches!(key, Key::Other(_)) {
        match raw_value {
            Some(MISSING) => None,
            Some(t) => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
            None => Some(Value::Flag),
        }
    } else if let Some(t) = raw_value {
        match t {
            MISSING => None,
            _ => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
        }
    } else {
        return Err(ParseError::MissingValue);
    };

    Ok((key, value))
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
