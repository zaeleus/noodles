pub(super) mod value;

use std::{error, fmt};

use self::value::parse_value;
use crate::{
    header::{record::value::map::info::Type, Number},
    io::reader::record::MISSING,
    record::info::field::Value,
    Header,
};

/// An error when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A value is missing.
    MissingValue(String),
    /// A value is invalid.
    InvalidValue(String, value::ParseError),
}

impl ParseError {
    /// Returns the key of the field that caused the failure.
    pub fn key(&self) -> Option<&str> {
        match self {
            Self::MissingValue(key) => Some(key),
            Self::InvalidValue(key, _) => Some(key),
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::MissingValue(_) => write!(f, "missing value"),
            ParseError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn parse_field(header: &Header, s: &str) -> Result<(String, Option<Value>), ParseError> {
    use crate::header::record::value::map::info::definition::definition;

    const MAX_COMPONENTS: usize = 2;
    const SEPARATOR: char = '=';

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let key = components.next().unwrap_or_default();

    let definition = header
        .infos()
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)));

    let raw_value = components.next();

    let value = if let Some((number, ty)) = definition {
        if matches!(ty, Type::Flag) {
            match raw_value.unwrap_or_default() {
                MISSING => None,
                t => parse_value(number, ty, t)
                    .map(Some)
                    .map_err(|e| ParseError::InvalidValue(key.into(), e))?,
            }
        } else if let Some(t) = raw_value {
            match t {
                MISSING => None,
                _ => parse_value(number, ty, t)
                    .map(Some)
                    .map_err(|e| ParseError::InvalidValue(key.into(), e))?,
            }
        } else {
            return Err(ParseError::MissingValue(key.into()));
        }
    } else {
        match raw_value {
            Some(MISSING) => None,
            Some(t) => parse_value(Number::Count(1), Type::String, t)
                .map(Some)
                .map_err(|e| ParseError::InvalidValue(key.into(), e))?,
            None => Some(Value::Flag),
        }
    };

    Ok((key.into(), value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::info::field::key;

    #[test]
    fn test_parse_field() {
        let header = Header::default();

        assert_eq!(
            parse_field(&header, "NS=2"),
            Ok((
                String::from(key::SAMPLES_WITH_DATA_COUNT),
                Some(Value::Integer(2))
            ))
        );

        assert!(matches!(
            parse_field(&header, "NS="),
            Err(ParseError::InvalidValue(s, _)) if s == key::SAMPLES_WITH_DATA_COUNT
        ));

        assert!(matches!(
            parse_field(&header, "NS=ndls"),
            Err(ParseError::InvalidValue(s, _)) if s == key::SAMPLES_WITH_DATA_COUNT
        ));
    }
}
