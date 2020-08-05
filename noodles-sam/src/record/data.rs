//! SAM record data and fields.

pub mod field;

pub use self::field::Field;

use std::{error, fmt, ops::Deref, str::FromStr};

const DELIMITER: char = '\t';

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data(Vec<Field>);

impl Deref for Data {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.iter().enumerate() {
            if i > 0 {
                f.write_str("\t")?;
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

impl From<Vec<Field>> for Data {
    fn from(fields: Vec<Field>) -> Self {
        Self(fields)
    }
}

/// An error returned when raw SAM record data fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// The input data contains an invalid field.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        s.split(DELIMITER)
            .map(|t| t.parse().map_err(ParseError::InvalidField))
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_fmt() {
        let data = Data::from(vec![
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
            Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
        ]);

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let data: Data = "RG:Z:rg0\tNH:i:1".parse()?;
        assert_eq!(data.len(), 2);

        assert!("".parse::<Data>()?.is_empty());

        Ok(())
    }
}
