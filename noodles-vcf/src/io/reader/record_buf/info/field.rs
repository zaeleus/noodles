pub(super) mod value;

use std::{error, fmt};

use self::value::parse_value;
use crate::{
    Header,
    header::record::value::map::info::{Number, Type},
    io::reader::record_buf::MISSING,
    variant::record_buf::info::field::Value,
};

/// An error when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The key is empty.
    EmptyKey,
    /// A value is missing.
    MissingValue(String),
    /// A value is invalid.
    InvalidValue(String, value::ParseError),
}

impl ParseError {
    /// Returns the key of the field that caused the failure.
    pub fn key(&self) -> Option<&str> {
        match self {
            Self::EmptyKey => None,
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
            ParseError::EmptyKey => write!(f, "empty key"),
            ParseError::MissingValue(_) => write!(f, "missing value"),
            ParseError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn parse_field(header: &Header, s: &str) -> Result<(String, Option<Value>), ParseError> {
    use crate::header::record::value::map::info::definition::definition;

    let (key, raw_value) = split_field(s)?;

    let definition = header
        .infos()
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)));

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

fn split_field(s: &str) -> Result<(&str, Option<&str>), ParseError> {
    const SEPARATOR: char = '=';

    let (key, value) = s
        .split_once(SEPARATOR)
        .map(|(k, v)| (k, Some(v)))
        .unwrap_or((s, None));

    if key.is_empty() {
        return Err(ParseError::EmptyKey);
    }

    Ok((key, value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::info::field::key;

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

    #[test]
    fn test_split_field() -> Result<(), ParseError> {
        assert_eq!(split_field("NS=2"), Ok(("NS", Some("2"))));
        assert_eq!(split_field("DP=."), Ok(("DP", Some("."))));
        assert_eq!(split_field("DP="), Ok(("DP", Some(""))));
        assert_eq!(split_field("H3"), Ok(("H3", None)));

        assert_eq!(split_field(""), Err(ParseError::EmptyKey));
        assert_eq!(split_field("="), Err(ParseError::EmptyKey));
        assert_eq!(split_field("=2"), Err(ParseError::EmptyKey));

        Ok(())
    }
}
