//! VCF record information and field.

pub mod field;

pub use self::field::Field;

use std::{error, fmt, str::FromStr};

use indexmap::IndexMap;

use super::MISSING_FIELD;
use crate::header;

const DELIMITER: char = ';';

/// VCF record information fields (`INFO`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Info(IndexMap<field::Key, Field>);

impl Info {
    /// Parses raw VCF record info.
    pub fn try_from_str(s: &str, infos: &header::Infos) -> Result<Self, ParseError> {
        parse(s, infos)
    }

    /// Returns the number of info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::Info;
    /// let info = Info::default();
    /// assert_eq!(info.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether there are any info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::Info;
    /// let info = Info::default();
    /// assert!(info.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns a reference to the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns, dp.clone()])?;
    ///
    /// assert_eq!(info.get(&Key::TotalDepth), Some(&dp));
    /// assert!(info.get(&Key::AlleleFrequencies).is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get(&self, key: &field::Key) -> Option<&Field> {
        self.0.get(key)
    }

    /// Returns a mutable reference to the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let mut info = Info::try_from(vec![ns, dp])?;
    ///
    /// if let Some(dp) = info.get_mut(&Key::TotalDepth) {
    ///     *dp.value_mut() = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(
    ///     info.get(&Key::TotalDepth).map(|field| field.value()),
    ///     Some(Some(&Value::Integer(8)))
    /// );
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_mut(&mut self, key: &field::Key) -> Option<&mut Field> {
        self.0.get_mut(key)
    }

    /// Returns a reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns, dp.clone()])?;
    ///
    /// assert_eq!(info.get_index(1), Some(&dp));
    /// assert!(info.get_index(5).is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_index(&self, i: usize) -> Option<&Field> {
        self.0.get_index(i).map(|(_, field)| field)
    }

    /// Returns a mutable reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let mut info = Info::try_from(vec![ns, dp])?;
    ///
    /// if let Some(dp) = info.get_index_mut(1) {
    ///     *dp.value_mut() = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(
    ///     info.get_index(1).map(|field| field.value()),
    ///     Some(Some(&Value::Integer(8)))
    /// );
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_index_mut(&mut self, i: usize) -> Option<&mut Field> {
        self.0.get_index_mut(i).map(|(_, field)| field)
    }

    /// Inserts a field into the info.
    ///
    /// This uses the field key as the key and field as the value.
    ///
    /// If the key already exists in the map, the existing field is replaced by the new one, and
    /// the existing field is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let mut info = Info::try_from(vec![ns])?;
    /// assert_eq!(info.len(), 1);
    ///
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// info.insert(dp.clone());
    ///
    /// assert_eq!(info.len(), 2);
    /// assert_eq!(info.get(&Key::TotalDepth), Some(&dp));
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn insert(&mut self, field: Field) -> Option<Field> {
        self.0.insert(field.key().clone(), field)
    }

    /// Returns an iterator over all keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns, dp.clone()])?;
    ///
    /// let mut keys = info.keys();
    ///
    /// assert_eq!(keys.next(), Some(&Key::SamplesWithDataCount));
    /// assert_eq!(keys.next(), Some(&Key::TotalDepth));
    /// assert!(keys.next().is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = &field::Key> {
        self.values().map(|field| field.key())
    }

    /// Returns an iterator over all fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{info::{field::{Key, Value}, Field}, Info};
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns.clone(), dp.clone()])?;
    ///
    /// let mut values = info.values();
    ///
    /// assert_eq!(values.next(), Some(&ns));
    /// assert_eq!(values.next(), Some(&dp));
    /// assert!(values.next().is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &Field> {
        self.0.values()
    }
}

impl AsRef<IndexMap<field::Key, Field>> for Info {
    fn as_ref(&self) -> &IndexMap<field::Key, Field> {
        &self.0
    }
}

impl AsMut<IndexMap<field::Key, Field>> for Info {
    fn as_mut(&mut self) -> &mut IndexMap<field::Key, Field> {
        &mut self.0
    }
}

impl fmt::Display for Info {
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
        Self::try_from_str(s, &header::Infos::default())
    }
}

fn parse(s: &str, infos: &header::Infos) -> Result<Info, ParseError> {
    match s {
        "" => Err(ParseError::Empty),
        MISSING_FIELD => Ok(Info::default()),
        _ => {
            let fields = s
                .split(DELIMITER)
                .map(|s| Field::try_from_str(s, infos))
                .collect::<Result<Vec<_>, _>>()
                .map_err(ParseError::InvalidField)?;

            Info::try_from(fields).map_err(ParseError::Invalid)
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
            Some(field::Value::Integer(2)),
        )])?;
        assert_eq!(info.to_string(), "NS=2");

        let info = Info::try_from(vec![
            Field::new(
                field::Key::SamplesWithDataCount,
                Some(field::Value::Integer(2)),
            ),
            Field::new(
                field::Key::AlleleFrequencies,
                Some(field::Value::FloatArray(vec![Some(0.333), Some(0.667)])),
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
            Some(field::Value::Integer(2)),
        )];
        assert!(Info::try_from(fields).is_ok());

        let fields = vec![
            Field::new(
                field::Key::SamplesWithDataCount,
                Some(field::Value::Integer(2)),
            ),
            Field::new(
                field::Key::SamplesWithDataCount,
                Some(field::Value::Integer(2)),
            ),
        ];
        assert_eq!(
            Info::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(
                field::Key::SamplesWithDataCount
            ))
        );
    }
}
