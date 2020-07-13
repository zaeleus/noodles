//! VCF record information and field.

pub mod field;

pub use self::field::Field;

use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ';';

/// VCF record information fields (`INFO`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Info(Vec<Field>);

impl Deref for Info {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, field) in self.iter().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?
                }

                write!(f, "{}", field)?;
            }

            Ok(())
        }
    }
}

impl From<Vec<Field>> for Info {
    fn from(fields: Vec<Field>) -> Self {
        Self(fields)
    }
}

/// An error returned when a raw VCF information fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
        }
    }
}

impl FromStr for Info {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Info::default()),
            _ => s
                .split(DELIMITER)
                .map(|s| s.parse())
                .collect::<Result<_, _>>()
                .map(Info)
                .map_err(ParseError::InvalidField),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let info = Info::default();
        assert_eq!(info.to_string(), ".");

        let info = Info(vec![Field::new(
            field::Key::SamplesWithDataCount,
            field::Value::Integer(2),
        )]);
        assert_eq!(info.to_string(), "NS=2");

        let info = Info(vec![
            Field::new(field::Key::SamplesWithDataCount, field::Value::Integer(2)),
            Field::new(
                field::Key::AlleleFrequencies,
                field::Value::FloatArray(vec![0.333, 0.667]),
            ),
        ]);
        assert_eq!(info.to_string(), "NS=2;AF=0.333,0.667");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Info = ".".parse()?;
        assert!(actual.is_empty());

        let actual: Info = "NS=2".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Info = "NS=2;AF=0.333,0.667".parse()?;
        assert_eq!(actual.len(), 2);

        assert_eq!("".parse::<Info>(), Err(ParseError::Empty));
        assert!(matches!(
            "NS=ndls".parse::<Info>(),
            Err(ParseError::InvalidField(_))
        ));

        Ok(())
    }
}
