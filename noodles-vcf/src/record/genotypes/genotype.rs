//! VCF record genotype and field.

pub mod field;

pub use self::field::Field;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
};

use indexmap::IndexMap;

use super::Keys;
use crate::{
    header::{
        format::Key,
        record::value::{map::Format, Map},
        Formats,
    },
    record::MISSING_FIELD,
};

const DELIMITER: char = ':';

/// A VCF record genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotype(IndexMap<Key, Field>);

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

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::Invalid(e) => Some(e),
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
            Self::InvalidField(_) => f.write_str("invalid field"),
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

impl error::Error for GenotypeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            Self::InvalidValueType(_) => None,
        }
    }
}

impl fmt::Display for GenotypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(_) => f.write_str("invalid value"),
            Self::InvalidValueType(value) => write!(f, "invalid String, got {value:?}"),
        }
    }
}

impl Genotype {
    /// Parses a raw genotype for the given genotype keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Key, record::value::{map::Format, Map}},
    ///     record::genotypes::{genotype::field::Value, Genotype},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Key::Genotype, Map::<Format>::from(&Key::Genotype))
    ///     .add_format(
    ///         Key::ConditionalGenotypeQuality,
    ///         Map::<Format>::from(&Key::ConditionalGenotypeQuality),
    ///     )
    ///     .build();
    ///
    /// let keys = "GT:GQ".parse()?;
    ///
    /// assert_eq!(
    ///     Genotype::parse("0|0:13", header.formats(), &keys),
    ///     Ok([
    ///         (Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         (Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ].into_iter().collect())
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn parse(s: &str, formats: &Formats, keys: &Keys) -> Result<Self, ParseError> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        } else if s == MISSING_FIELD {
            return Ok(Self::default());
        }

        let mut fields = Vec::with_capacity(keys.len());

        for (raw_field, key) in s.split(DELIMITER).zip(keys.iter()) {
            let field = if let Some(format) = formats.get(key) {
                Field::from_str_format(raw_field, key, format).map_err(ParseError::InvalidField)?
            } else {
                let format = Map::<Format>::from(key);
                Field::from_str_format(raw_field, key, &format).map_err(ParseError::InvalidField)?
            };

            fields.push(field);
        }

        Self::try_from(fields).map_err(ParseError::Invalid)
    }

    /// Returns the VCF record genotypes genotype value.
    ///
    /// This is a convenience method to return a parsed version of the genotype (`GT`) field value.
    /// Use `[Self::get]` with `[Key::Genotype]` to get the raw value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Key, record::value::{map::Format, Map}},
    ///     record::genotypes::{genotype::{field::Value, Field}, Genotype},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Key::Genotype, Map::<Format>::from(&Key::Genotype))
    ///     .add_format(
    ///         Key::ConditionalGenotypeQuality,
    ///         Map::<Format>::from(&Key::ConditionalGenotypeQuality),
    ///     )
    ///     .build();
    ///
    /// let keys = "GT:GQ".parse()?;
    ///
    /// let genotype = Genotype::parse("0|0:13", header.formats(), &keys)?;
    /// assert_eq!(genotype.genotype(), Some(Ok("0|0".parse()?)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotype(&self) -> Option<Result<field::value::Genotype, GenotypeError>> {
        self.get(&Key::Genotype)
            .and_then(|f| f.value())
            .map(|value| match value {
                field::Value::String(s) => s.parse().map_err(GenotypeError::InvalidValue),
                _ => Err(GenotypeError::InvalidValueType(value.clone())),
            })
    }
}

impl Deref for Genotype {
    type Target = IndexMap<Key, Field>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genotype {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, field) in self.values().enumerate() {
                if i > 0 {
                    write!(f, "{DELIMITER}")?;
                }

                write!(f, "{field}")?;
            }

            Ok(())
        }
    }
}

impl Extend<(Key, Option<field::Value>)> for Genotype {
    fn extend<T: IntoIterator<Item = (Key, Option<field::Value>)>>(&mut self, iter: T) {
        for (key, value) in iter {
            let field = Field::new(key.clone(), value);
            self.insert(key, field);
        }
    }
}

impl FromIterator<(Key, Option<field::Value>)> for Genotype {
    fn from_iter<T: IntoIterator<Item = (Key, Option<field::Value>)>>(iter: T) -> Self {
        let mut genotype = Self::default();
        genotype.extend(iter);
        genotype
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
    DuplicateKey(Key),
}

impl error::Error for TryFromFieldsError {}

impl fmt::Display for TryFromFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGenotypeFieldPosition => f.write_str("invalid genotype field position"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

impl TryFrom<Vec<Field>> for Genotype {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        if fields.is_empty() {
            return Ok(Self::default());
        } else if let Some(i) = fields.iter().position(|f| f.key() == &Key::Genotype) {
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
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::Header::builder()
            .add_format(Key::Genotype, Map::<Format>::from(&Key::Genotype))
            .add_format(
                Key::ConditionalGenotypeQuality,
                Map::<Format>::from(&Key::ConditionalGenotypeQuality),
            )
            .build();

        let keys = "GT".parse()?;
        assert_eq!(
            Genotype::parse(".", header.formats(), &keys),
            Ok(Genotype::default())
        );

        let keys = "GT".parse()?;
        let actual = Genotype::parse("0|0", header.formats(), &keys)?;
        assert_eq!(actual.len(), 1);

        let keys = "GT:GQ".parse()?;
        let actual = Genotype::parse("0|0:13", header.formats(), &keys)?;
        assert_eq!(actual.len(), 2);

        let keys = "GT".parse()?;
        assert_eq!(
            Genotype::parse("", header.formats(), &keys),
            Err(ParseError::Empty)
        );

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let genotype = Genotype::default();
        assert_eq!(genotype.to_string(), ".");

        let genotype: Genotype = [(
            Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )]
        .into_iter()
        .collect();

        assert_eq!(genotype.to_string(), "0|0");

        let genotype: Genotype = [
            (
                Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
            (
                Key::ConditionalGenotypeQuality,
                Some(field::Value::Integer(13)),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(genotype.to_string(), "0|0:13");
    }

    #[test]
    fn test_extend() {
        let mut genotype = Genotype::default();

        let fields = [(
            Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )];
        genotype.extend(fields);

        let expected = [(
            Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )]
        .into_iter()
        .collect();

        assert_eq!(genotype, expected);
    }

    #[test]
    fn test_try_from_fields_for_genotype() {
        assert!(Genotype::try_from(Vec::new()).is_ok());

        let fields = vec![Field::new(
            Key::Genotype,
            Some(field::Value::String(String::from("0|0"))),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![Field::new(
            Key::ConditionalGenotypeQuality,
            Some(field::Value::Integer(13)),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![
            Field::new(
                Key::ConditionalGenotypeQuality,
                Some(field::Value::Integer(13)),
            ),
            Field::new(
                Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::InvalidGenotypeFieldPosition)
        );

        let fields = vec![
            Field::new(
                Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
            Field::new(
                Key::Genotype,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(Key::Genotype))
        );
    }

    #[test]
    fn test_genotype() {
        let genotype: Genotype = [(
            Key::Genotype,
            Some(field::Value::String(String::from("ndls"))),
        )]
        .into_iter()
        .collect();

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValue(_)))
        ));

        let genotype: Genotype = [(Key::Genotype, Some(field::Value::Integer(0)))]
            .into_iter()
            .collect();

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValueType(_)))
        ));
    }
}
