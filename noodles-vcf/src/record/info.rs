//! VCF record information and field.

pub mod field;

use std::{error, fmt, hash::Hash, str::FromStr};

use indexmap::IndexMap;

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

    /// Removes all fields from the info map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let mut info: Info = [ns, dp].into_iter().collect();
    /// assert!(!info.is_empty());
    ///
    /// info.clear();
    /// assert!(info.is_empty());
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let info: Info = [ns, dp.clone()].into_iter().collect();
    ///
    /// assert_eq!(info.get(&key::TOTAL_DEPTH), Some(Some(&Value::Integer(13))));
    /// assert!(info.get(&key::ALLELE_FREQUENCIES).is_none());
    /// ```
    pub fn get<K>(&self, key: &K) -> Option<Option<&field::Value>>
    where
        K: Hash + indexmap::Equivalent<Key>,
    {
        self.0.get(key).map(|value| value.as_ref())
    }

    /// Returns a mutable reference to the field value with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let mut info: Info = [ns, dp].into_iter().collect();
    ///
    /// if let Some(value) = info.get_mut(&key::TOTAL_DEPTH) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(info.get(&key::TOTAL_DEPTH), Some(Some(&Value::Integer(8))));
    /// ```
    pub fn get_mut<K>(&mut self, key: &K) -> Option<&mut Option<field::Value>>
    where
        K: Hash + indexmap::Equivalent<Key>,
    {
        self.0.get_mut(key)
    }

    /// Returns a reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// assert_eq!(
    ///     info.get_index(1),
    ///     Some((&key::TOTAL_DEPTH, Some(&Value::Integer(13))))
    /// );
    ///
    /// assert!(info.get_index(5).is_none());
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let mut info: Info = [ns, dp].into_iter().collect();
    ///
    /// if let Some((_, value)) = info.get_index_mut(1) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(
    ///     info.get_index(1),
    ///     Some((&key::TOTAL_DEPTH, Some(&Value::Integer(8))))
    /// );
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let mut info: Info = [ns].into_iter().collect();
    /// assert_eq!(info.len(), 1);
    ///
    /// info.insert(key::TOTAL_DEPTH, Some(Value::Integer(13)));
    ///
    /// assert_eq!(info.len(), 2);
    /// assert_eq!(info.get(&key::TOTAL_DEPTH), Some(Some(&Value::Integer(13))));
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// let mut keys = info.keys();
    ///
    /// assert_eq!(keys.next(), Some(&key::SAMPLES_WITH_DATA_COUNT));
    /// assert_eq!(keys.next(), Some(&key::TOTAL_DEPTH));
    /// assert!(keys.next().is_none());
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = &Key> {
        self.0.keys()
    }

    /// Returns an iterator over all values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::info::key,
    ///     record::{info::field::Value, Info},
    /// };
    ///
    /// let ns = (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)));
    /// let dp = (key::TOTAL_DEPTH, Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// let mut values = info.values();
    ///
    /// assert_eq!(values.next(), Some(Some(&Value::Integer(2))));
    /// assert_eq!(values.next(), Some(Some(&Value::Integer(13))));
    /// assert!(values.next().is_none());
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
        for (i, (key, value)) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, "{DELIMITER}")?;
            }

            key.fmt(f)?;

            match value {
                None => f.write_str("=.")?,
                Some(field::Value::Flag) => {}
                Some(v) => write!(f, "={v}")?,
            }
        }

        Ok(())
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
        self.0.extend(iter);
    }
}

impl FromIterator<(Key, Option<field::Value>)> for Info {
    fn from_iter<T: IntoIterator<Item = (Key, Option<field::Value>)>>(iter: T) -> Self {
        let mut info = Self::default();
        info.extend(iter);
        info
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
        _ => {
            let mut info = Info::default();

            for raw_field in s.split(DELIMITER) {
                let (key, value) =
                    field::parse(raw_field, infos).map_err(ParseError::InvalidField)?;

                if info.insert(key.clone(), value).is_some() {
                    return Err(ParseError::Invalid(TryFromFieldsError::DuplicateKey(key)));
                }
            }

            Ok(info)
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
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::info::key;

    #[test]
    fn test_fmt() {
        let info = Info::default();
        assert!(info.to_string().is_empty());

        let info: Info = [(key::SAMPLES_WITH_DATA_COUNT, Some(field::Value::Integer(2)))]
            .into_iter()
            .collect();
        assert_eq!(info.to_string(), "NS=2");

        let info: Info = [
            (key::SAMPLES_WITH_DATA_COUNT, Some(field::Value::Integer(2))),
            (
                key::ALLELE_FREQUENCIES,
                Some(field::Value::FloatArray(vec![Some(0.333), Some(0.667)])),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(info.to_string(), "NS=2;AF=0.333,0.667");
    }

    #[test]
    fn test_extend() {
        let mut info = Info::default();

        let fields = [(key::SAMPLES_WITH_DATA_COUNT, Some(field::Value::Integer(2)))];
        info.extend(fields);

        let expected = [(key::SAMPLES_WITH_DATA_COUNT, Some(field::Value::Integer(2)))]
            .into_iter()
            .collect();

        assert_eq!(info, expected);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Info = "NS=2".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Info = "NS=2;AF=0.333,0.667".parse()?;
        assert_eq!(actual.len(), 2);

        assert_eq!("".parse::<Info>(), Err(ParseError::Empty));
        assert!(matches!(
            ".".parse::<Info>(),
            Err(ParseError::InvalidField(_))
        ));
        assert!(matches!(
            "NS=ndls".parse::<Info>(),
            Err(ParseError::InvalidField(_))
        ));

        Ok(())
    }
}
