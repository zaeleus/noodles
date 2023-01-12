//! VCF record information and field.

pub mod field;

pub use self::field::Field;

use std::{error, fmt, str::FromStr};

use indexmap::IndexMap;

use super::MISSING_FIELD;
use crate::header::{self, info::Key};

const DELIMITER: char = ';';

/// VCF record information fields (`INFO`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Info(IndexMap<Key, Option<field::Value>>);

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

    /// Removes all field from the info map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let mut info = Info::try_from(vec![ns, dp])?;
    /// assert!(!info.is_empty());
    ///
    /// info.clear();
    /// assert!(info.is_empty());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Returns a reference to the field value with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns, dp.clone()])?;
    ///
    /// assert_eq!(info.get(&Key::TotalDepth), Some(dp.value()));
    /// assert!(info.get(&Key::AlleleFrequencies).is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get(&self, key: &Key) -> Option<Option<&field::Value>> {
        self.0.get(key).map(|value| value.as_ref())
    }

    /// Returns a mutable reference to the field value with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let mut info = Info::try_from(vec![ns, dp])?;
    ///
    /// if let Some(value) = info.get_mut(&Key::TotalDepth) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(info.get(&Key::TotalDepth), Some(Some(&Value::Integer(8))));
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_mut(&mut self, key: &Key) -> Option<&mut Option<field::Value>> {
        self.0.get_mut(key)
    }

    /// Returns a reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns, dp.clone()])?;
    ///
    /// assert_eq!(info.get_index(1), Some((dp.key(), dp.value())));
    /// assert!(info.get_index(5).is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_index(&self, i: usize) -> Option<(&Key, Option<&field::Value>)> {
        self.0
            .get_index(i)
            .map(|(key, value)| (key, value.as_ref()))
    }

    /// Returns a mutable reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let mut info = Info::try_from(vec![ns, dp])?;
    ///
    /// if let Some((_, value)) = info.get_index_mut(1) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(
    ///     info.get_index(1),
    ///     Some((&Key::TotalDepth, Some(&Value::Integer(8))))
    /// );
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn get_index_mut(&mut self, i: usize) -> Option<(&mut Key, &mut Option<field::Value>)> {
        self.0.get_index_mut(i)
    }

    /// Inserts a field into the info map.
    ///
    /// If the key already exists in the map, the existing value is replaced by the new one, and
    /// the existing value is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let mut info = Info::try_from(vec![ns])?;
    /// assert_eq!(info.len(), 1);
    ///
    /// info.insert(Key::TotalDepth, Some(Value::Integer(13)));
    ///
    /// assert_eq!(info.len(), 2);
    /// assert_eq!(info.get(&Key::TotalDepth), Some(Some(&Value::Integer(13))));
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn insert(
        &mut self,
        key: Key,
        value: Option<field::Value>,
    ) -> Option<Option<field::Value>> {
        self.0.insert(key, value)
    }

    /// Returns an iterator over all keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
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
    pub fn keys(&self) -> impl Iterator<Item = &Key> {
        self.0.keys()
    }

    /// Returns an iterator over all fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::Key,
    ///     record::{info::{field::Value, Field}, Info},
    /// };
    ///
    /// let ns = Field::new(Key::SamplesWithDataCount, Some(Value::Integer(2)));
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// let info = Info::try_from(vec![ns.clone(), dp.clone()])?;
    ///
    /// let mut values = info.values();
    ///
    /// assert_eq!(values.next(), Some(ns.value()));
    /// assert_eq!(values.next(), Some(dp.value()));
    /// assert!(values.next().is_none());
    /// # Ok::<_, noodles_vcf::record::info::TryFromFieldsError>(())
    /// ```
    pub fn values(&self) -> impl Iterator<Item = Option<&field::Value>> {
        self.0.values().map(|value| value.as_ref())
    }
}

impl AsRef<IndexMap<Key, Option<field::Value>>> for Info {
    fn as_ref(&self) -> &IndexMap<Key, Option<field::Value>> {
        &self.0
    }
}

impl AsMut<IndexMap<Key, Option<field::Value>>> for Info {
    fn as_mut(&mut self) -> &mut IndexMap<Key, Option<field::Value>> {
        &mut self.0
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, (key, value)) in self.0.iter().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?;
                }

                key.fmt(f)?;

                match value {
                    None => f.write_str("=.")?,
                    Some(field::Value::Flag) => {}
                    Some(v) => write!(f, "={}", v)?,
                }
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

impl Extend<(Key, Option<field::Value>)> for Info {
    fn extend<T: IntoIterator<Item = (Key, Option<field::Value>)>>(&mut self, iter: T) {
        for (key, value) in iter {
            self.insert(key, value);
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
    DuplicateKey(Key),
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
            let (key, value) = (field.key().clone(), field.value().cloned());

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
    fn test_fmt() -> Result<(), TryFromFieldsError> {
        let info = Info::default();
        assert_eq!(info.to_string(), ".");

        let info = Info::try_from(vec![Field::new(
            Key::SamplesWithDataCount,
            Some(field::Value::Integer(2)),
        )])?;
        assert_eq!(info.to_string(), "NS=2");

        let info = Info::try_from(vec![
            Field::new(Key::SamplesWithDataCount, Some(field::Value::Integer(2))),
            Field::new(
                Key::AlleleFrequencies,
                Some(field::Value::FloatArray(vec![Some(0.333), Some(0.667)])),
            ),
        ])?;
        assert_eq!(info.to_string(), "NS=2;AF=0.333,0.667");

        Ok(())
    }

    #[test]
    fn test_extend() -> Result<(), TryFromFieldsError> {
        let mut info = Info::default();

        let fields = [(Key::SamplesWithDataCount, Some(field::Value::Integer(2)))];
        info.extend(fields);

        let expected = Info::try_from(vec![Field::new(
            Key::SamplesWithDataCount,
            Some(field::Value::Integer(2)),
        )])?;

        assert_eq!(info, expected);

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
            Key::SamplesWithDataCount,
            Some(field::Value::Integer(2)),
        )];
        assert!(Info::try_from(fields).is_ok());

        let fields = vec![
            Field::new(Key::SamplesWithDataCount, Some(field::Value::Integer(2))),
            Field::new(Key::SamplesWithDataCount, Some(field::Value::Integer(2))),
        ];
        assert_eq!(
            Info::try_from(fields),
            Err(TryFromFieldsError::DuplicateKey(Key::SamplesWithDataCount))
        );
    }
}
