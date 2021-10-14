//! VCF record genotype and field.

pub mod field;

pub use self::field::Field;

use std::{convert::TryFrom, error, fmt, ops::Deref};

use indexmap::IndexMap;

use crate::record::{Format, MISSING_FIELD};

const DELIMITER: char = ':';

/// A VCF record genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotype(IndexMap<field::Key, Field>);

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
            Self::Invalid(e) => write!(f, "invalid input: {}", e),
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
        }
    }
}

/// An error returned when a genotype (`GT`) field value fails to parse.
#[derive(Clone, Debug, PartialEq)]
pub enum GenotypeError {
    /// The genotype field value is invalid.
    InvalidValue(field::value::genotype::ParseError),
    /// The genotype field value type is invalid.
    ///
    /// The `GT` field value must be a `String`.
    InvalidValueType(field::Value),
}

impl error::Error for GenotypeError {}

impl fmt::Display for GenotypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
            Self::InvalidValueType(value) => write!(f, "invalid String, got {:?}", value),
        }
    }
}

impl Genotype {
    /// Parses a raw genotype for the given genotype format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::record::genotypes::{genotype::{field::{Key, Value}, Field}, Genotype};
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

    /// Returns the VCF record genotypes genotype value.
    ///
    /// This is a convenience method to return a parsed version of the genotype (`GT`) field value.
    /// Use `[Self::get]` with `[field::Key::Genotype]` to get the raw value.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::record::{genotypes::{genotype::{field::{Key, Value}, Field}, Genotype}};
    /// let genotype = Genotype::from_str_format("0|0:13", &"GT:GQ".parse()?)?;
    /// assert_eq!(genotype.genotype(), Some(Ok("0|0".parse()?)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotype(&self) -> Option<Result<field::value::Genotype, GenotypeError>> {
        self.get(&field::Key::Genotype)
            .and_then(|f| f.value())
            .map(|value| match value {
                field::Value::String(s) => s.parse().map_err(GenotypeError::InvalidValue),
                _ => Err(GenotypeError::InvalidValueType(value.clone())),
            })
    }
}

impl Deref for Genotype {
    type Target = IndexMap<field::Key, Field>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, field) in self.values().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?;
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
    /// A key is duplicated.
    ///
    /// § 1.6.2 Genotype fields (2021-01-13): "...duplicate keys are not allowed".
    DuplicateKey(field::Key),
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGenotypeFieldPosition => f.write_str("invalid genotype field position"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {}", key),
        }
    }
}

impl TryFrom<Vec<Field>> for Genotype {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        if fields.is_empty() {
            return Ok(Self::default());
        } else if let Some(i) = fields.iter().position(|f| f.key() == &field::Key::Genotype) {
            if i != 0 {
                return Err(TryFromFieldsError::InvalidGenotypeFieldPosition);
            }
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
    fn test_from_str_format() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT".parse()?;
        assert_eq!(
            Genotype::from_str_format(".", &format),
            Ok(Genotype::default())
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
    fn test_fmt() -> Result<(), TryFromFieldsError> {
        let genotype = Genotype::default();
        assert_eq!(genotype.to_string(), ".");

        let genotype = Genotype::try_from(vec![Field::new(
            field::Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )])?;

        assert_eq!(genotype.to_string(), "0|0");

        let genotype = Genotype::try_from(vec![
            Field::new(
                field::Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
            Field::new(
                field::Key::ConditionalGenotypeQuality,
                Some(field::Value::Integer(13)),
            ),
        ])?;

        assert_eq!(genotype.to_string(), "0|0:13");

        Ok(())
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

        let fields = vec![
            Field::new(
                field::Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
            Field::new(
                field::Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(field::Key::Genotype))
        );
    }

    #[test]
    fn test_genotype() -> Result<(), TryFromFieldsError> {
        let genotype = Genotype::try_from(vec![Field::new(
            field::Key::Genotype,
            Some(field::Value::String(String::from("ndls"))),
        )])?;

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValue(_)))
        ));

        let genotype = Genotype::try_from(vec![Field::new(
            field::Key::Genotype,
            Some(field::Value::Integer(0)),
        )])?;

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValueType(_)))
        ));

        Ok(())
    }
}
