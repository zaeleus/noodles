//! VCF record information and field.

pub mod field;

use std::hash::Hash;

use indexmap::IndexMap;

use self::field::Value;

/// VCF record information fields (`INFO`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Info(IndexMap<String, Option<Value>>);

impl Info {
    /// Removes all fields from the info map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::Info as _,
    ///     record_buf::{info::field::{key, Value}, Info },
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
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
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let info: Info = [ns, dp.clone()].into_iter().collect();
    ///
    /// assert_eq!(info.get(key::TOTAL_DEPTH), Some(Some(&Value::Integer(13))));
    /// assert!(info.get(key::ALLELE_FREQUENCIES).is_none());
    /// ```
    pub fn get<K>(&self, key: &K) -> Option<Option<&Value>>
    where
        K: Hash + indexmap::Equivalent<String> + ?Sized,
    {
        self.0.get(key).map(|value| value.as_ref())
    }

    /// Returns a mutable reference to the field value with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let mut info: Info = [ns, dp].into_iter().collect();
    ///
    /// if let Some(value) = info.get_mut(key::TOTAL_DEPTH) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(info.get(key::TOTAL_DEPTH), Some(Some(&Value::Integer(8))));
    /// ```
    pub fn get_mut<K>(&mut self, key: &K) -> Option<&mut Option<Value>>
    where
        K: Hash + indexmap::Equivalent<String> + ?Sized,
    {
        self.0.get_mut(key)
    }

    /// Returns a reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// assert_eq!(
    ///     info.get_index(1),
    ///     Some((&String::from(key::TOTAL_DEPTH), Some(&Value::Integer(13))))
    /// );
    ///
    /// assert!(info.get_index(5).is_none());
    /// ```
    pub fn get_index(&self, i: usize) -> Option<(&String, Option<&Value>)> {
        self.0
            .get_index(i)
            .map(|(key, value)| (key, value.as_ref()))
    }

    /// Returns a mutable reference to the field at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let mut info: Info = [ns, dp].into_iter().collect();
    ///
    /// if let Some((_, value)) = info.get_index_mut(1) {
    ///     *value = Some(Value::Integer(8));
    /// }
    ///
    /// assert_eq!(
    ///     info.get_index(1),
    ///     Some((&String::from(key::TOTAL_DEPTH), Some(&Value::Integer(8))))
    /// );
    /// ```
    pub fn get_index_mut(&mut self, i: usize) -> Option<(&String, &mut Option<Value>)> {
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
    /// use noodles_vcf::variant::{
    ///     record::Info as _,
    ///     record_buf::{info::field::{key, Value}, Info},
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let mut info: Info = [ns].into_iter().collect();
    /// assert_eq!(info.len(), 1);
    ///
    /// info.insert(String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    ///
    /// assert_eq!(info.len(), 2);
    /// assert_eq!(info.get(key::TOTAL_DEPTH), Some(Some(&Value::Integer(13))));
    /// ```
    pub fn insert(&mut self, key: String, value: Option<Value>) -> Option<Option<Value>> {
        self.0.insert(key, value)
    }

    /// Returns an iterator over all keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// let mut keys = info.keys();
    ///
    /// assert_eq!(keys.next(), Some(&String::from(key::SAMPLES_WITH_DATA_COUNT)));
    /// assert_eq!(keys.next(), Some(&String::from(key::TOTAL_DEPTH)));
    /// assert!(keys.next().is_none());
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = &String> {
        self.0.keys()
    }

    /// Returns an iterator over all values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{
    ///     info::field::{key, Value},
    ///     Info,
    /// };
    ///
    /// let ns = (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(2)));
    /// let dp = (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    /// let info: Info = [ns, dp].into_iter().collect();
    ///
    /// let mut values = info.values();
    ///
    /// assert_eq!(values.next(), Some(Some(&Value::Integer(2))));
    /// assert_eq!(values.next(), Some(Some(&Value::Integer(13))));
    /// assert!(values.next().is_none());
    /// ```
    pub fn values(&self) -> impl Iterator<Item = Option<&Value>> {
        self.0.values().map(|value| value.as_ref())
    }
}

impl AsRef<IndexMap<String, Option<Value>>> for Info {
    fn as_ref(&self) -> &IndexMap<String, Option<Value>> {
        &self.0
    }
}

impl AsMut<IndexMap<String, Option<Value>>> for Info {
    fn as_mut(&mut self) -> &mut IndexMap<String, Option<Value>> {
        &mut self.0
    }
}

impl Extend<(String, Option<Value>)> for Info {
    fn extend<T: IntoIterator<Item = (String, Option<Value>)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(String, Option<Value>)> for Info {
    fn from_iter<T: IntoIterator<Item = (String, Option<Value>)>>(iter: T) -> Self {
        let mut info = Self::default();
        info.extend(iter);
        info
    }
}

impl crate::variant::record::Info for Info {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        _: &'h crate::Header,
    ) -> Box<
        dyn Iterator<
                Item = std::io::Result<(
                    &'a str,
                    Option<crate::variant::record::info::field::Value<'a>>,
                )>,
            > + 'a,
    > {
        Box::new(
            self.0
                .iter()
                .map(|(key, value)| Ok((key.as_ref(), value.as_ref().map(|v| v.into())))),
        )
    }
}

impl crate::variant::record::Info for &Info {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        _: &'h crate::Header,
    ) -> Box<
        dyn Iterator<
                Item = std::io::Result<(
                    &'a str,
                    Option<crate::variant::record::info::field::Value<'a>>,
                )>,
            > + 'a,
    > {
        Box::new(
            self.0
                .iter()
                .map(|(key, value)| Ok((key.as_ref(), value.as_ref().map(|v| v.into())))),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{field::key, *};

    #[test]
    fn test_extend() {
        let mut info = Info::default();

        let fields = [(
            String::from(key::SAMPLES_WITH_DATA_COUNT),
            Some(Value::from(2)),
        )];
        info.extend(fields);

        let expected = [(
            String::from(key::SAMPLES_WITH_DATA_COUNT),
            Some(Value::from(2)),
        )]
        .into_iter()
        .collect();

        assert_eq!(info, expected);
    }
}
