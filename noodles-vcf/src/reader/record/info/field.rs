pub(super) mod value;

use std::{error, fmt};

use self::value::parse_value;
use crate::{
    header::record::value::map::info::Type,
    reader::record::MISSING,
    record::info::field::{key, Key, Value},
    Header,
};

/// An error when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A key is invalid.
    InvalidKey(key::ParseError),
    /// A value is missing.
    MissingValue(Key),
    /// A value is invalid.
    InvalidValue(Key, value::ParseError),
}

impl ParseError {
    /// Returns the key of the field that caused the failure.
    pub fn key(&self) -> Option<&Key> {
        match self {
            Self::InvalidKey(_) => None,
            Self::MissingValue(key) => Some(key),
            Self::InvalidValue(key, _) => Some(key),
        }
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidKey(_) => write!(f, "invalid key"),
            ParseError::MissingValue(_) => write!(f, "missing value"),
            ParseError::InvalidValue(..) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn parse_field(header: &Header, s: &str) -> Result<(Key, Option<Value>), ParseError> {
    use crate::header::record::value::map::info::definition::definition;

    const MAX_COMPONENTS: usize = 2;
    const SEPARATOR: char = '=';

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let raw_key = components.next().unwrap_or_default();
    let key = raw_key.parse().map_err(ParseError::InvalidKey)?;

    let (number, ty) = header
        .infos()
        .get(&key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), &key).map(|(n, t, _)| (n, t)))
        .unwrap_or_default();

    let raw_value = components.next();

    let value = if matches!(ty, Type::Flag) {
        match raw_value.unwrap_or_default() {
            MISSING => None,
            t => parse_value(number, ty, t)
                .map(Some)
                .map_err(|e| ParseError::InvalidValue(key.clone(), e))?,
        }
    } else if matches!(key, Key::Other(_)) {
        match raw_value {
            Some(MISSING) => None,
            Some(t) => parse_value(number, ty, t)
                .map(Some)
                .map_err(|e| ParseError::InvalidValue(key.clone(), e))?,
            None => Some(Value::Flag),
        }
    } else if let Some(t) = raw_value {
        match t {
            MISSING => None,
            _ => parse_value(number, ty, t)
                .map(Some)
                .map_err(|e| ParseError::InvalidValue(key.clone(), e))?,
        }
    } else {
        return Err(ParseError::MissingValue(key));
    };

    Ok((key, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() {
        let header = Header::default();

        assert_eq!(
            parse_field(&header, "NS=2"),
            Ok((key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2))))
        );

        assert!(matches!(
            parse_field(&header, "."),
            Err(ParseError::InvalidKey(_))
        ));

        assert!(matches!(
            parse_field(&header, "NS="),
            Err(ParseError::InvalidValue(key::SAMPLES_WITH_DATA_COUNT, _))
        ));

        assert!(matches!(
            parse_field(&header, "NS=ndls"),
            Err(ParseError::InvalidValue(key::SAMPLES_WITH_DATA_COUNT, _))
        ));
    }
}
