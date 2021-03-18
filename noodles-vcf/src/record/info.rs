//! VCF record information and field.

pub mod field;

pub use self::field::Field;

use std::{convert::TryFrom, error, fmt, ops::Deref, str::FromStr};

use indexmap::IndexMap;

use super::MISSING_FIELD;

const DELIMITER: char = ';';

/// VCF record information fields (`INFO`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Info(IndexMap<field::Key, Field>);

impl Deref for Info {
    type Target = IndexMap<field::Key, Field>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, field) in self.values().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?
                }

                write!(f, "{}", field)?;
            }

            Ok(())
        }
    }
}

/// An error returned when a raw VCF information fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(TryFromFieldsError),
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(e) => write!(f, "invalid input: {}", e),
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
        }
    }
}

impl FromStr for Info {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self::default()),
            _ => {
                let fields = s
                    .split(DELIMITER)
                    .map(|s| s.parse())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(ParseError::InvalidField)?;

                Self::try_from(fields).map_err(ParseError::Invalid)
            }
        }
    }
}

/// An error returned when VCF info fields fail to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFieldsError {
    /// A key is duplicated.
    ///
    /// § 1.6.1 Fixed fields (2021-01-13): "Duplicate keys are not allowed."
    DuplicateKey(field::Key),
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateKey(key) => write!(f, "duplicate key: {}", key),
        }
    }
}

impl TryFrom<Vec<Field>> for Info {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        if fields.is_empty() {
            return Ok(Self::default());
        }

        let mut map = IndexMap::with_capacity(fields.len());

        for field in fields {
            if let Some(duplicate_field) = map.insert(field.key().clone(), field) {
                return Err(TryFromFieldsError::DuplicateKey(
                    duplicate_field.key().clone(),
                ));
            }
        }

        Ok(Self(map))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), TryFromFieldsError> {
        let info = Info::default();
        assert_eq!(info.to_string(), ".");

        let info = Info::try_from(vec![Field::new(
            field::Key::SamplesWithDataCount,
            field::Value::Integer(2),
        )])?;
        assert_eq!(info.to_string(), "NS=2");

        let info = Info::try_from(vec![
            Field::new(field::Key::SamplesWithDataCount, field::Value::Integer(2)),
            Field::new(
                field::Key::AlleleFrequencies,
                field::Value::FloatArray(vec![0.333, 0.667]),
            ),
        ])?;
        assert_eq!(info.to_string(), "NS=2;AF=0.333,0.667");

        Ok(())
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

    #[test]
    fn test_try_from_fields_for_info() {
        assert_eq!(Info::try_from(Vec::new()), Ok(Info::default()));

        let fields = vec![Field::new(
            field::Key::SamplesWithDataCount,
            field::Value::Integer(2),
        )];
        assert!(Info::try_from(fields).is_ok());

        let fields = vec![
            Field::new(field::Key::SamplesWithDataCount, field::Value::Integer(2)),
            Field::new(field::Key::SamplesWithDataCount, field::Value::Integer(2)),
        ];
        assert_eq!(
            Info::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(
                field::Key::SamplesWithDataCount
            ))
        );
    }
}
