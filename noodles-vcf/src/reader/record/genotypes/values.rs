use std::{error, fmt};

use noodles_core as core;

use crate::{
    header::record::value::{map::Format, Map},
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
    InvalidValue(crate::record::genotypes::sample::value::ParseError),
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
    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == MISSING {
        return Ok(());
    }

    let mut raw_values = s.split(DELIMITER);

    for (key, raw_value) in keys.iter().zip(&mut raw_values) {
        let value = if let Some(format) = header.formats().get(key) {
            parse_value(format, raw_value).map_err(ParseError::InvalidValue)?
        } else {
            let format = Map::<Format>::from(key);
            parse_value(&format, raw_value).map_err(ParseError::InvalidValue)?
        };

        values.push(value);
    }

    if raw_values.next().is_some() {
        return Err(ParseError::UnexpectedValue);
    }

    Ok(())
}

fn parse_value(
    format: &Map<Format>,
    s: &str,
) -> Result<Option<Value>, crate::record::genotypes::sample::value::ParseError> {
    match s {
        MISSING => Ok(None),
        _ => Value::from_str_format(s, format).map(Some),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_values() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::format::key;

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
