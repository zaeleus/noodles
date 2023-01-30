//! VCF record info field.

pub mod value;

pub use self::value::Value;

use std::{error, fmt};

use crate::header::{
    self,
    info::{Key, Type},
    record::value::{map::Info, Map},
    Infos,
};

const MISSING_VALUE: &str = ".";
const SEPARATOR: char = '=';

/// An error returned when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The key is missing.
    MissingKey,
    /// The key is invalid.
    InvalidKey(header::info::key::ParseError),
    /// The value is missing.
    MissingValue,
    /// The value is invalid.
    InvalidValue(value::ParseError),
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
            Self::MissingKey => f.write_str("missing key"),
            Self::InvalidKey(_) => f.write_str("invalid key"),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
        }
    }
}

pub(super) fn parse(s: &str, infos: &Infos) -> Result<(Key, Option<Value>), ParseError> {
    const MAX_COMPONENTS: usize = 2;

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let key: Key = components
        .next()
        .ok_or(ParseError::MissingKey)
        .and_then(|t| t.parse().map_err(ParseError::InvalidKey))?;

    let value = if let Some(info) = infos.get(&key) {
        parse_value(&mut components, &key, info)?
    } else {
        let info = Map::<Info>::from(&key);
        parse_value(&mut components, &key, &info)?
    };

    Ok((key, value))
}

fn parse_value<'a, I>(
    iter: &mut I,
    key: &Key,
    info: &Map<Info>,
) -> Result<Option<Value>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    if let Type::Flag = info.ty() {
        let t = iter.next().unwrap_or_default();

        if t == MISSING_VALUE {
            Ok(None)
        } else {
            Value::from_str_info(t, info)
                .map(Some)
                .map_err(ParseError::InvalidValue)
        }
    } else if matches!(key, Key::Other(..)) {
        if let Some(t) = iter.next() {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                Value::from_str_info(t, info)
                    .map(Some)
                    .map_err(ParseError::InvalidValue)
            }
        } else {
            Ok(Some(Value::Flag))
        }
    } else if let Some(t) = iter.next() {
        if t == MISSING_VALUE {
            Ok(None)
        } else {
            Value::from_str_info(t, info)
                .map(Some)
                .map_err(ParseError::InvalidValue)
        }
    } else {
        Err(ParseError::MissingValue)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), crate::header::info::key::ParseError> {
        let header = crate::Header::builder()
            .add_info(Key::AlleleCount, Map::<Info>::from(&Key::AlleleCount))
            .add_info(
                Key::SamplesWithDataCount,
                Map::<Info>::from(&Key::SamplesWithDataCount),
            )
            .add_info(
                Key::IsSomaticMutation,
                Map::<Info>::from(&Key::IsSomaticMutation),
            )
            .add_info(
                Key::BreakendEventId,
                Map::<Info>::from(&Key::BreakendEventId),
            )
            .build();

        assert_eq!(parse("AC=.", header.infos()), Ok((Key::AlleleCount, None)));

        assert_eq!(
            parse("NS=2", header.infos()),
            Ok((Key::SamplesWithDataCount, Some(Value::Integer(2))))
        );

        assert_eq!(
            parse("BQ=1.333", header.infos()),
            Ok((Key::BaseQuality, Some(Value::Float(1.333))))
        );

        assert_eq!(
            parse("SOMATIC", header.infos()),
            Ok((Key::IsSomaticMutation, Some(Value::Flag)))
        );

        assert_eq!(
            parse("EVENT=INV0", header.infos()),
            Ok((
                Key::BreakendEventId,
                Some(Value::String(String::from("INV0")))
            ))
        );

        let key = "NDLS".parse()?;
        assert_eq!(
            parse("NDLS=VCF", header.infos()),
            Ok((key, Some(Value::String(String::from("VCF")))))
        );

        let key = "FLG".parse()?;
        assert_eq!(parse("FLG", header.infos()), Ok((key, Some(Value::Flag))));

        Ok(())
    }
}
