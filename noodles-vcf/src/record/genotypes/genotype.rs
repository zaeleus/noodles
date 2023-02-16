//! VCF record genotype and field.

pub mod field;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
};

use indexmap::IndexMap;

use super::Keys;
use crate::{
    header::{
        format::{key, Key},
        record::value::{map::Format, Map},
        Formats,
    },
    record::MISSING_FIELD,
};

const DELIMITER: char = ':';

/// A VCF record genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotype(IndexMap<Key, Option<field::Value>>);

/// An error returned when a raw VCF genotype fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(TryFromFieldsError),
    /// A field is invalid.
    InvalidField(field::ParseError),
    /// The genotype field value is unexpected.
    ///
    /// There are fewer keys than values.
    UnexpectedValue,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            Self::InvalidField(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
            Self::InvalidField(_) => f.write_str("invalid field"),
            Self::UnexpectedValue => f.write_str("unexpected value"),
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
    InvalidValueType(Option<field::Value>),
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
    ///     header::{format::key, record::value::{map::Format, Map}},
    ///     record::genotypes::{genotype::field::Value, Genotype},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
    ///     .add_format(
    ///         key::CONDITIONAL_GENOTYPE_QUALITY,
    ///         Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
    ///     )
    ///     .build();
    ///
    /// let keys = "GT:GQ".parse()?;
    ///
    /// assert_eq!(
    ///     Genotype::parse("0|0:13", header.formats(), &keys),
    ///     Ok([
    ///         (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
    ///         (key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13))),
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
        let mut raw_values = s.split(DELIMITER);

        for (key, raw_value) in keys.iter().zip(&mut raw_values) {
            let field = if let Some(format) = formats.get(key) {
                let value = parse_value(format, raw_value).map_err(ParseError::InvalidField)?;
                (key.clone(), value)
            } else {
                let format = Map::<Format>::from(key);
                let value = parse_value(&format, raw_value).map_err(ParseError::InvalidField)?;
                (key.clone(), value)
            };

            fields.push(field);
        }

        if raw_values.next().is_some() {
            Err(ParseError::UnexpectedValue)
        } else {
            Self::try_from(fields).map_err(ParseError::Invalid)
        }
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
    ///     header::{format::key, record::value::{map::Format, Map}},
    ///     record::genotypes::{genotype::field::Value, Genotype},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
    ///     .add_format(
    ///         key::CONDITIONAL_GENOTYPE_QUALITY,
    ///         Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
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
        self.get(&key::GENOTYPE).map(|value| match value {
            Some(field::Value::String(s)) => s.parse().map_err(GenotypeError::InvalidValue),
            _ => Err(GenotypeError::InvalidValueType(value.clone())),
        })
    }
}

impl Deref for Genotype {
    type Target = IndexMap<Key, Option<field::Value>>;

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
            for (i, value) in self.values().enumerate() {
                if i > 0 {
                    write!(f, "{DELIMITER}")?;
                }

                if let Some(v) = value {
                    write!(f, "{v}")?;
                } else {
                    f.write_str(".")?;
                }
            }

            Ok(())
        }
    }
}

impl Extend<(Key, Option<field::Value>)> for Genotype {
    fn extend<T: IntoIterator<Item = (Key, Option<field::Value>)>>(&mut self, iter: T) {
        self.0.extend(iter);
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

impl TryFrom<Vec<(Key, Option<field::Value>)>> for Genotype {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<(Key, Option<field::Value>)>) -> Result<Self, Self::Error> {
        if fields.is_empty() {
            return Ok(Self::default());
        } else if let Some(i) = fields.iter().position(|(key, _)| key == &key::GENOTYPE) {
            if i != 0 {
                return Err(TryFromFieldsError::InvalidGenotypeFieldPosition);
            }
        }

        let mut map = IndexMap::with_capacity(fields.len());

        for (key, value) in fields {
            if map.insert(key.clone(), value).is_some() {
                return Err(TryFromFieldsError::DuplicateKey(key));
            }
        }

        Ok(Self(map))
    }
}

fn parse_value(format: &Map<Format>, s: &str) -> Result<Option<field::Value>, field::ParseError> {
    if s == "." {
        Ok(None)
    } else {
        field::Value::from_str_format(s, format)
            .map(Some)
            .map_err(field::ParseError::InvalidValue)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::Header::builder()
            .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
            .add_format(
                key::CONDITIONAL_GENOTYPE_QUALITY,
                Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
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

        let keys = "GT".parse()?;
        assert_eq!(
            Genotype::parse("0|0:13", header.formats(), &keys),
            Err(ParseError::UnexpectedValue)
        );

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let genotype = Genotype::default();
        assert_eq!(genotype.to_string(), ".");

        let genotype: Genotype = [(
            key::GENOTYPE,
            Some(field::Value::String(String::from("0|0"))),
        )]
        .into_iter()
        .collect();

        assert_eq!(genotype.to_string(), "0|0");

        let genotype: Genotype = [
            (
                key::GENOTYPE,
                Some(field::Value::String(String::from("0|0"))),
            ),
            (
                key::CONDITIONAL_GENOTYPE_QUALITY,
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
            key::GENOTYPE,
            Some(field::Value::String(String::from("0|0"))),
        )];
        genotype.extend(fields);

        let expected = [(
            key::GENOTYPE,
            Some(field::Value::String(String::from("0|0"))),
        )]
        .into_iter()
        .collect();

        assert_eq!(genotype, expected);
    }

    #[test]
    fn test_try_from_fields_for_genotype() {
        assert!(Genotype::try_from(Vec::new()).is_ok());

        let fields = vec![(
            key::GENOTYPE,
            Some(field::Value::String(String::from("0|0"))),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![(
            key::CONDITIONAL_GENOTYPE_QUALITY,
            Some(field::Value::Integer(13)),
        )];
        assert!(Genotype::try_from(fields).is_ok());

        let fields = vec![
            (
                key::CONDITIONAL_GENOTYPE_QUALITY,
                Some(field::Value::Integer(13)),
            ),
            (
                key::GENOTYPE,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::InvalidGenotypeFieldPosition)
        );

        let fields = vec![
            (
                key::GENOTYPE,
                Some(field::Value::String(String::from("0|0"))),
            ),
            (
                key::GENOTYPE,
                Some(field::Value::String(String::from("0|0"))),
            ),
        ];
        assert_eq!(
            Genotype::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(key::GENOTYPE))
        );
    }

    #[test]
    fn test_genotype() {
        let genotype: Genotype = [(
            key::GENOTYPE,
            Some(field::Value::String(String::from("ndls"))),
        )]
        .into_iter()
        .collect();

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValue(_)))
        ));

        let genotype: Genotype = [(key::GENOTYPE, Some(field::Value::Integer(0)))]
            .into_iter()
            .collect();

        assert!(matches!(
            genotype.genotype(),
            Some(Err(GenotypeError::InvalidValueType(_)))
        ));
    }
}
