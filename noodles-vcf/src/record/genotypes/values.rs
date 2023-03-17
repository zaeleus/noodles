//! VCF record genotype values.

pub mod value;

pub use self::value::Value;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
};

use indexmap::IndexMap;

use super::Keys;
use crate::{
    header::{
        format::{key, Key},
        Formats,
    },
    record::MISSING_FIELD,
};

const DELIMITER: char = ':';

/// A VCF record genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Values(IndexMap<Key, Option<Value>>);

/// An error returned when a raw VCF genotype fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(TryFromFieldsError),
    /// A value is invalid.
    InvalidValue(value::ParseError),
    /// The genotype field value is unexpected.
    ///
    /// There are fewer keys than values.
    UnexpectedValue,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
            Self::UnexpectedValue => f.write_str("unexpected value"),
        }
    }
}

/// An error returned when a genotype (`GT`) field value fails to parse.
#[derive(Clone, Debug, PartialEq)]
pub enum GenotypeError {
    /// The genotype field value is invalid.
    InvalidValue(value::genotype::ParseError),
    /// The genotype field value type is invalid.
    ///
    /// The `GT` field value must be a `String`.
    InvalidValueType(Option<Value>),
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

impl Values {
    /// Parses a raw genotype for the given genotype keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::key, record::value::{map::Format, Map}},
    ///     record::genotypes::{values::Value, Values},
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
    ///     Values::parse("0|0:13", header.formats(), &keys),
    ///     Ok([
    ///         (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
    ///         (key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13))),
    ///     ].into_iter().collect())
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn parse(s: &str, formats: &Formats, keys: &Keys) -> Result<Self, ParseError> {
        super::parse_values(s, formats, keys)
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
    ///     record::genotypes::{values::Value, Values},
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
    /// let values = Values::parse("0|0:13", header.formats(), &keys)?;
    /// assert_eq!(values.genotype(), Some(Ok("0|0".parse()?)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotype(&self) -> Option<Result<value::Genotype, GenotypeError>> {
        self.get(&key::GENOTYPE).map(|value| match value {
            Some(Value::String(s)) => s.parse().map_err(GenotypeError::InvalidValue),
            _ => Err(GenotypeError::InvalidValueType(value.clone())),
        })
    }
}

impl Deref for Values {
    type Target = IndexMap<Key, Option<Value>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Values {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Values {
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

impl Extend<(Key, Option<Value>)> for Values {
    fn extend<T: IntoIterator<Item = (Key, Option<Value>)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(Key, Option<Value>)> for Values {
    fn from_iter<T: IntoIterator<Item = (Key, Option<Value>)>>(iter: T) -> Self {
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

impl TryFrom<Vec<(Key, Option<Value>)>> for Values {
    type Error = TryFromFieldsError;

    fn try_from(fields: Vec<(Key, Option<Value>)>) -> Result<Self, Self::Error> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::record::value::{map::Format, Map};

        let header = crate::Header::builder()
            .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
            .add_format(
                key::CONDITIONAL_GENOTYPE_QUALITY,
                Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
            )
            .build();

        let keys = "GT".parse()?;
        assert_eq!(
            Values::parse(".", header.formats(), &keys),
            Ok(Values::default())
        );

        let keys = "GT".parse()?;
        let actual = Values::parse("0|0", header.formats(), &keys)?;
        assert_eq!(actual.len(), 1);

        let keys = "GT:GQ".parse()?;
        let actual = Values::parse("0|0:13", header.formats(), &keys)?;
        assert_eq!(actual.len(), 2);

        let keys = "GT".parse()?;
        assert_eq!(
            Values::parse("", header.formats(), &keys),
            Err(ParseError::Empty)
        );

        let keys = "GT".parse()?;
        assert_eq!(
            Values::parse("0|0:13", header.formats(), &keys),
            Err(ParseError::UnexpectedValue)
        );

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let values = Values::default();
        assert_eq!(values.to_string(), ".");

        let values: Values = [(key::GENOTYPE, Some(Value::String(String::from("0|0"))))]
            .into_iter()
            .collect();

        assert_eq!(values.to_string(), "0|0");

        let values: Values = [
            (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
            (key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13))),
        ]
        .into_iter()
        .collect();

        assert_eq!(values.to_string(), "0|0:13");
    }

    #[test]
    fn test_extend() {
        let mut values = Values::default();

        let fields = [(key::GENOTYPE, Some(Value::String(String::from("0|0"))))];
        values.extend(fields);

        let expected = [(key::GENOTYPE, Some(Value::String(String::from("0|0"))))]
            .into_iter()
            .collect();

        assert_eq!(values, expected);
    }

    #[test]
    fn test_try_from_fields_for_values() {
        assert!(Values::try_from(Vec::new()).is_ok());

        let fields = vec![(key::GENOTYPE, Some(Value::String(String::from("0|0"))))];
        assert!(Values::try_from(fields).is_ok());

        let fields = vec![(key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13)))];
        assert!(Values::try_from(fields).is_ok());

        let fields = vec![
            (key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13))),
            (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
        ];
        assert_eq!(
            Values::try_from(fields),
            Err(TryFromFieldsError::InvalidGenotypeFieldPosition)
        );

        let fields = vec![
            (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
            (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
        ];
        assert_eq!(
            Values::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(key::GENOTYPE))
        );
    }

    #[test]
    fn test_genotype() {
        let values: Values = [(key::GENOTYPE, Some(Value::String(String::from("ndls"))))]
            .into_iter()
            .collect();

        assert!(matches!(
            values.genotype(),
            Some(Err(GenotypeError::InvalidValue(_)))
        ));

        let values: Values = [(key::GENOTYPE, Some(Value::Integer(0)))]
            .into_iter()
            .collect();

        assert!(matches!(
            values.genotype(),
            Some(Err(GenotypeError::InvalidValueType(_)))
        ));
    }
}
