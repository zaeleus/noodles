//! VCF record info field.

pub mod key;
pub mod value;

pub use self::{key::Key, value::Value};

use std::{error, fmt};

use crate::header::{
    record::value::{
        map::{info::Type, Info},
        Map,
    },
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
    InvalidKey(key::ParseError),
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

pub(crate) fn parse(s: &str, infos: &Infos) -> Result<(Key, Option<Value>), ParseError> {
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

pub(crate) fn parse_value<'a, I>(
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
    fn test_parse() -> Result<(), key::ParseError> {
        let header = crate::Header::builder()
            .add_info(key::ALLELE_COUNT, Map::<Info>::from(&key::ALLELE_COUNT))
            .add_info(
                key::SAMPLES_WITH_DATA_COUNT,
                Map::<Info>::from(&key::SAMPLES_WITH_DATA_COUNT),
            )
            .add_info(
                key::IS_SOMATIC_MUTATION,
                Map::<Info>::from(&key::IS_SOMATIC_MUTATION),
            )
            .add_info(
                key::BREAKEND_EVENT_ID,
                Map::<Info>::from(&key::BREAKEND_EVENT_ID),
            )
            .build();

        assert_eq!(parse("AC=.", header.infos()), Ok((key::ALLELE_COUNT, None)));

        assert_eq!(
            parse("NS=2", header.infos()),
            Ok((key::SAMPLES_WITH_DATA_COUNT, Some(Value::from(2))))
        );

        assert_eq!(
            parse("BQ=1.333", header.infos()),
            Ok((key::BASE_QUALITY, Some(Value::from(1.333))))
        );

        assert_eq!(
            parse("SOMATIC", header.infos()),
            Ok((key::IS_SOMATIC_MUTATION, Some(Value::Flag)))
        );

        assert_eq!(
            parse("EVENT=INV0", header.infos()),
            Ok((
                key::BREAKEND_EVENT_ID,
                Some(Value::from(vec![Some(String::from("INV0"))]))
            ))
        );

        let key = "NDLS".parse()?;
        assert_eq!(
            parse("NDLS=VCF", header.infos()),
            Ok((key, Some(Value::from("VCF"))))
        );

        let key = "FLG".parse()?;
        assert_eq!(parse("FLG", header.infos()), Ok((key, Some(Value::Flag))));

        Ok(())
    }
}
