//! VCF record genotype and field.

pub mod field;

pub use self::field::Field;

use std::{convert::TryFrom, error, fmt, ops::Deref};

use super::{Format, MISSING_FIELD};

const DELIMITER: char = ':';

/// A VCF record genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotype(Vec<Field>);

/// An error returned when a raw VCF genotype fails to parse.
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
            Self::Invalid(e) => write!(f, "{}", e),
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
        }
    }
}

impl Genotype {
    /// Parses a raw genotype for the given genotype format.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::convert::TryFrom;
    /// use noodles_vcf::record::{genotype::{field::{Key, Value}, Field}, Genotype};
    ///
    /// let format = "GT:GQ".parse()?;
    ///
    /// assert_eq!(
    ///     Genotype::from_str_format("0|0:13", &format),
    ///     Ok(Genotype::try_from(vec![
    ///         Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ])?)
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_str_format(s: &str, format: &Format) -> Result<Self, ParseError> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self::default()),
            _ => {
                let fields = s
                    .split(DELIMITER)
                    .zip(format.iter())
                    .map(|(t, k)| Field::from_str_key(t, k))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(ParseError::InvalidField)?;

                Self::try_from(fields).map_err(ParseError::Invalid)
            }
        }
    }
}

impl Deref for Genotype {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Genotype {
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

/// An error returned when VCF genotype fields fail to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromFieldsError {
    /// The genotype field (`GT`) position is invalid.
    ///
    /// The genotype field must be first if present. See § 1.6.2 Genotype fields (2020-06-25).
    InvalidGenotypeFieldPosition,
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGenotypeFieldPosition => f.write_str("invalid genotype field position"),
        }
    }
}

impl TryFrom<Vec<Field>> for Genotype {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        if fields.is_empty() {
            Ok(Self::default())
        } else if let Some(i) = fields.iter().position(|f| f.key() == &field::Key::Genotype) {
            if i == 0 {
                Ok(Self(fields))
            } else {
                Err(TryFromFieldsError::InvalidGenotypeFieldPosition)
            }
        } else {
            Ok(Self(fields))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_format() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT".parse()?;
        assert_eq!(
            Genotype::from_str_format(".", &format),
            Ok(Genotype(Vec::new()))
        );

        let format = "GT".parse()?;
        let actual = Genotype::from_str_format("0|0", &format)?;
        assert_eq!(actual.len(), 1);

        let format = "GT:GQ".parse()?;
        let actual = Genotype::from_str_format("0|0:13", &format)?;
        assert_eq!(actual.len(), 2);

        let format = "GT".parse()?;
        assert_eq!(
            Genotype::from_str_format("", &format),
            Err(ParseError::Empty)
        );

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let genotype = Genotype::default();
        assert_eq!(genotype.to_string(), ".");

        let genotype = Genotype(vec![Field::new(
            field::Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )]);

        assert_eq!(genotype.to_string(), "0|0");

        let genotype = Genotype(vec![
            Field::new(
                field::Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
            Field::new(
                field::Key::ConditionalGenotypeQuality,
                Some(field::Value::Integer(13)),
            ),
        ]);

        assert_eq!(genotype.to_string(), "0|0:13");
    }

    #[test]
    fn test_try_from_fields_for_genotype() {
        assert!(Genotype::try_from(Vec::new()).is_ok());

        let fields = vec![Field::new(
            field::Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![Field::new(
            field::Key::ConditionalGenotypeQuality,
            Some(field::Value::Integer(13)),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![
            Field::new(
                field::Key::ConditionalGenotypeQuality,
                Some(field::Value::Integer(13)),
            ),
            Field::new(
                field::Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::InvalidGenotypeFieldPosition)
        );
    }
}
