//! Alignment record data buffer.

pub mod field;

use std::{io, mem};

use self::field::Value;
use crate::alignment::record::data::field::Tag;

/// An alignment record data buffer.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Data(Vec<(Tag, Value)>);

impl Data {
    /// Returns the number of fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::Data;
    /// let data = Data::default();
    /// assert_eq!(data.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether there are any fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::Data;
    /// let data = Data::default();
    /// assert!(data.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Removes all fields from the data map.
    ///
    /// This does not affect the internal capacity.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let nh = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let mut data: Data = [nh].into_iter().collect();
    /// assert_eq!(data.len(), 1);
    /// data.clear();
    /// assert!(data.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Returns a reference to the value of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let (tag, value) = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// assert_eq!(data.get(&tag), Some(&value));
    /// assert!(data.get(&Tag::READ_GROUP).is_none());
    /// ```
    pub fn get<K>(&self, tag: &K) -> Option<&Value>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.0
            .iter()
            .find(|(t, _)| tag.equivalent(t))
            .map(|(_, v)| v)
    }

    /// Returns the index of the field of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let nh = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [nh].into_iter().collect();
    ///
    /// assert_eq!(data.get_index_of(&Tag::ALIGNMENT_HIT_COUNT), Some(0));
    /// assert!(data.get_index_of(&Tag::READ_GROUP).is_none());
    /// ```
    pub fn get_index_of<K>(&self, tag: &K) -> Option<usize>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.0.iter().position(|(t, _)| tag.equivalent(t))
    }

    /// Returns an iterator over all tag-value pairs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let (tag, value) = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// let mut fields = data.iter();
    /// assert_eq!(fields.next(), Some((tag, &value)));
    /// assert!(fields.next().is_none());
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (Tag, &Value)> {
        self.0.iter().map(|(tag, value)| (*tag, value))
    }

    /// Returns an iterator over all tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let nh = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [nh].into_iter().collect();
    ///
    /// let mut keys = data.keys();
    /// assert_eq!(keys.next(), Some(Tag::ALIGNMENT_HIT_COUNT));
    /// assert!(keys.next().is_none());
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = Tag> + '_ {
        self.0.iter().map(|(tag, _)| *tag)
    }

    /// Returns an iterator over all values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let (tag, value) = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// let mut values = data.values();
    /// assert_eq!(values.next(), Some(&value));
    /// assert!(values.next().is_none());
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &Value> {
        self.0.iter().map(|(_, value)| value)
    }

    /// Inserts a field into the data map.
    ///
    /// This uses the field tag as the key and field as the value.
    ///
    /// If the tag already exists in the map, the existing field is replaced by the new one, and
    /// the existing field is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let mut data = Data::default();
    /// data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: Value) -> Option<(Tag, Value)> {
        let field = (tag, value);

        match self.get_index_of(&tag) {
            Some(i) => Some(mem::replace(&mut self.0[i], field)),
            None => {
                self.0.push(field);
                None
            }
        }
    }

    /// Removes the field with the given tag.
    ///
    /// The field is returned if it exists.
    ///
    /// This works like [`Vec::swap_remove`]; it does not preserve the order but has a constant
    /// time complexity.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::Tag,
    ///     record_buf::{data::field::Value, Data},
    /// };
    ///
    /// let nh = (Tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let rg = (Tag::READ_GROUP, Value::from("rg0"));
    /// let md = (Tag::ALIGNMENT_SCORE, Value::from(98));
    /// let mut data: Data = [nh.clone(), rg.clone(), md.clone()].into_iter().collect();
    ///
    /// assert_eq!(data.remove(&Tag::ALIGNMENT_HIT_COUNT), Some(nh));
    /// assert!(data.remove(&Tag::COMMENT).is_none());
    ///
    /// let expected = [md, rg].into_iter().collect();
    /// assert_eq!(data, expected);
    /// ```
    pub fn remove<K>(&mut self, tag: &K) -> Option<(Tag, Value)>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.swap_remove(tag)
    }

    fn swap_remove<K>(&mut self, tag: &K) -> Option<(Tag, Value)>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.get_index_of(tag).map(|i| self.0.swap_remove(i))
    }
}

impl crate::alignment::record::field::Data for &Data {
    fn is_empty(&self) -> bool {
        Data::is_empty(self)
    }

    fn get(
        &self,
        tag: &Tag,
    ) -> Option<io::Result<crate::alignment::record::data::field::Value<'_>>> {
        Data::get(self, tag).map(|value| Ok(value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<(Tag, crate::alignment::record::data::field::Value<'_>)>>
            + '_,
    > {
        Box::new(Data::iter(self).map(|(tag, value)| Ok((tag, value.into()))))
    }
}

impl crate::alignment::record::field::Data for Data {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &Tag,
    ) -> Option<io::Result<crate::alignment::record::data::field::Value<'_>>> {
        self.get(tag).map(|value| Ok(value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<(Tag, crate::alignment::record::data::field::Value<'_>)>>
            + '_,
    > {
        Box::new(self.iter().map(|(tag, value)| Ok((tag, value.into()))))
    }
}

impl Extend<(Tag, Value)> for Data {
    fn extend<T: IntoIterator<Item = (Tag, Value)>>(&mut self, iter: T) {
        for (tag, value) in iter {
            self.insert(tag, value);
        }
    }
}

impl FromIterator<(Tag, Value)> for Data {
    fn from_iter<T: IntoIterator<Item = (Tag, Value)>>(iter: T) -> Self {
        let mut data = Self::default();
        data.extend(iter);
        data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remove_with_multiple_removes() {
        let zz = Tag::new(b'z', b'z');

        let mut data: Data = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(2)),
            (Tag::EDIT_DISTANCE, Value::from(1)),
            (zz, Value::from(0)),
        ]
        .into_iter()
        .collect();

        data.remove(&Tag::EDIT_DISTANCE);
        data.remove(&zz);
        data.remove(&Tag::ALIGNMENT_HIT_COUNT);

        assert!(data.is_empty());
    }
}
