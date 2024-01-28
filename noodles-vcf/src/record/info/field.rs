//! VCF record info field.

pub mod key;
pub mod value;

pub use self::value::Value;

use std::{error, fmt};

use crate::header::{
    record::value::{
        map::{info::Type, Info},
        Map,
    },
    FileFormat, Infos, Number,
};

const MISSING_VALUE: &str = ".";
const SEPARATOR: char = '=';

/// An error returned when a raw VCF record info field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The key is missing.
    MissingKey,
    /// The value is missing.
    MissingValue,
    /// The value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingKey => f.write_str("missing key"),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
        }
    }
}

pub(crate) fn parse(s: &str, infos: &Infos) -> Result<(String, Option<Value>), ParseError> {
    const MAX_COMPONENTS: usize = 2;

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);
    let key = components.next().ok_or(ParseError::MissingKey)?;
    let value = parse_value(&mut components, infos, key)?;
    Ok((key.into(), value))
}

pub(crate) fn parse_value<'a, I>(
    iter: &mut I,
    infos: &Infos,
    key: &str,
) -> Result<Option<Value>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    use crate::header::record::value::map::info::definition::definition;

    let definition = infos
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(FileFormat::default(), key).map(|(n, t, _)| (n, t)));

    if let Some((number, ty)) = definition {
        let info = Map::<Info>::new(number, ty, "");

        if let Type::Flag = ty {
            let t = iter.next().unwrap_or_default();

            if t == MISSING_VALUE {
                Ok(None)
            } else {
                Value::from_str_info(t, &info)
                    .map(Some)
                    .map_err(ParseError::InvalidValue)
            }
        } else if let Some(t) = iter.next() {
            if t == MISSING_VALUE {
                Ok(None)
            } else {
                Value::from_str_info(t, &info)
                    .map(Some)
                    .map_err(ParseError::InvalidValue)
            }
        } else {
            Err(ParseError::MissingValue)
        }
    } else if let Some(t) = iter.next() {
        let info = Map::<Info>::new(Number::Count(1), Type::String, "");

        if t == MISSING_VALUE {
            Ok(None)
        } else {
            Value::from_str_info(t, &info)
                .map(Some)
                .map_err(ParseError::InvalidValue)
        }
    } else {
        Ok(Some(Value::Flag))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        let header = crate::Header::builder()
            .add_info(key::ALLELE_COUNT, Map::<Info>::from(key::ALLELE_COUNT))
            .add_info(
                key::SAMPLES_WITH_DATA_COUNT,
                Map::<Info>::from(key::SAMPLES_WITH_DATA_COUNT),
            )
            .add_info(
                key::IS_SOMATIC_MUTATION,
                Map::<Info>::from(key::IS_SOMATIC_MUTATION),
            )
            .add_info(
                key::BREAKEND_EVENT_ID,
                Map::<Info>::from(key::BREAKEND_EVENT_ID),
            )
            .build();

        assert_eq!(
            parse("AC=.", header.infos()),
            Ok((String::from(key::ALLELE_COUNT), None))
        );

        assert_eq!(
            parse("NS=2", header.infos()),
            Ok((
                String::from(key::SAMPLES_WITH_DATA_COUNT),
                Some(Value::from(2))
            ))
        );

        assert_eq!(
            parse("BQ=1.333", header.infos()),
            Ok((String::from(key::BASE_QUALITY), Some(Value::from(1.333))))
        );

        assert_eq!(
            parse("SOMATIC", header.infos()),
            Ok((String::from(key::IS_SOMATIC_MUTATION), Some(Value::Flag)))
        );

        assert_eq!(
            parse("EVENT=INV0", header.infos()),
            Ok((
                String::from(key::BREAKEND_EVENT_ID),
                Some(Value::from(vec![Some(String::from("INV0"))]))
            ))
        );

        assert_eq!(
            parse("NDLS=VCF", header.infos()),
            Ok((String::from("NDLS"), Some(Value::from("VCF"))))
        );

        assert_eq!(
            parse("FLG", header.infos()),
            Ok((String::from("FLG"), Some(Value::Flag)))
        );
    }
}
