mod value;

use std::{error, fmt};

use noodles_core as core;

use self::value::parse_value;
use crate::{
    header::{record::value::map::format::Type, Number},
    reader::record::MISSING,
    record::genotypes::{sample::Value, Keys},
    Header,
};

/// An error when raw VCF record genotypes values fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A value is invalid.
    InvalidValue(value::ParseError),
    /// The value was unexpected.
    ///
    /// There are unexpectedly more values than keys.
    UnexpectedValue,
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
            Self::Empty => f.write_str("empty input"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
            Self::UnexpectedValue => f.write_str("unexpected value"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_values(
    header: &Header,
    keys: &Keys,
    s: &str,
    values: &mut Vec<Option<Value>>,
) -> Result<(), ParseError> {
    use crate::header::record::value::map::format::definition::definition;

    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == MISSING {
        return Ok(());
    }

    let mut raw_values = s.split(DELIMITER);

    for (key, raw_value) in keys.iter().zip(&mut raw_values) {
        let value = match raw_value {
            MISSING => None,
            _ => {
                let (number, ty) = header
                    .formats()
                    .get(key)
                    .map(|format| (format.number(), format.ty()))
                    .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)))
                    .unwrap_or((Number::Count(1), Type::String));

                parse_value(number, ty, raw_value)
                    .map(Some)
                    .map_err(ParseError::InvalidValue)?
            }
        };

        values.push(value);
    }

    if raw_values.next().is_some() {
        return Err(ParseError::UnexpectedValue);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_values() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::genotypes::keys::key;

        let header = Header::default();
        let mut values = Vec::new();

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        parse_values(&header, &keys, ".", &mut values)?;
        assert!(values.is_empty());

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        parse_values(&header, &keys, "0|0", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0")))]);

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0:13", &mut values)?;
        assert_eq!(
            values,
            vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]
        );

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0:.", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0"))), None]);

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0")))]);

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        assert_eq!(
            parse_values(&header, &keys, "", &mut values),
            Err(ParseError::Empty)
        );

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        assert_eq!(
            parse_values(&header, &keys, "0|0:13", &mut values),
            Err(ParseError::UnexpectedValue)
        );

        Ok(())
    }
}
