//! SAM record data and fields.

pub mod field;

use std::{io, mem};

use self::field::{Tag, Value};

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Data {
    fields: Vec<(Tag, Value)>,
}

impl Data {
    /// Returns the number of data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Data;
    /// let data = Data::default();
    /// assert_eq!(data.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// Returns whether there are any data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Data;
    /// let data = Data::default();
    /// assert!(data.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// Removes all data fields from the data map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    /// let nh = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let mut data: Data = [nh].into_iter().collect();
    /// assert_eq!(data.len(), 1);
    /// data.clear();
    /// assert!(data.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.fields.clear();
    }

    /// Returns a reference to the value of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let (tag, value) = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// assert_eq!(data.get(&tag), Some(&value));
    /// assert!(data.get(&tag::READ_GROUP).is_none());
    /// ```
    pub fn get<K>(&self, tag: &K) -> Option<&Value>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.fields
            .iter()
            .find(|(t, _)| tag.equivalent(t))
            .map(|(_, v)| v)
    }

    /// Returns the index of the field of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let nh = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [nh].into_iter().collect();
    ///
    /// assert_eq!(data.get_index_of(&tag::ALIGNMENT_HIT_COUNT), Some(0));
    /// assert!(data.get_index_of(&tag::READ_GROUP).is_none());
    /// ```
    pub fn get_index_of<K>(&self, tag: &K) -> Option<usize>
    where
        K: indexmap::Equivalent<Tag>,
    {
        self.fields.iter().position(|(t, _)| tag.equivalent(t))
    }

    /// Returns an iterator over all tag-value pairs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let (tag, value) = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// let mut fields = data.iter();
    /// assert_eq!(fields.next(), Some((tag, &value)));
    /// assert!(fields.next().is_none());
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (Tag, &Value)> {
        self.fields.iter().map(|(tag, value)| (*tag, value))
    }

    /// Returns an iterator over all tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let nh = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [nh].into_iter().collect();
    ///
    /// let mut keys = data.keys();
    /// assert_eq!(keys.next(), Some(tag::ALIGNMENT_HIT_COUNT));
    /// assert!(keys.next().is_none());
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = Tag> + '_ {
        self.fields.iter().map(|(tag, _)| *tag)
    }

    /// Returns an iterator over all values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let (tag, value) = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let data: Data = [(tag, value.clone())].into_iter().collect();
    ///
    /// let mut values = data.values();
    /// assert_eq!(values.next(), Some(&value));
    /// assert!(values.next().is_none());
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &Value> {
        self.fields.iter().map(|(_, value)| value)
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
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    /// let mut data = Data::default();
    /// data.insert(tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: Value) -> Option<(Tag, Value)> {
        let field = (tag, value);

        match self.get_index_of(&tag) {
            Some(i) => Some(mem::replace(&mut self.fields[i], field)),
            None => {
                self.fields.push(field);
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
    /// use noodles_sam::record::{data::field::{tag, Value}, Data};
    ///
    /// let nh = (tag::ALIGNMENT_HIT_COUNT, Value::from(1));
    /// let rg = (tag::READ_GROUP, Value::from("rg0"));
    /// let md = (tag::ALIGNMENT_SCORE, Value::from(98));
    /// let mut data: Data = [nh.clone(), rg.clone(), md.clone()].into_iter().collect();
    ///
    /// assert_eq!(data.remove(&tag::ALIGNMENT_HIT_COUNT), Some(nh));
    /// assert!(data.remove(&tag::COMMENT).is_none());
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
        self.get_index_of(tag).map(|i| self.fields.swap_remove(i))
    }
}

impl crate::alignment::record::Data for &Data {
    fn is_empty(&self) -> bool {
        Data::is_empty(self)
    }

    fn get(
        &self,
        tag: &[u8; 2],
    ) -> Option<io::Result<crate::alignment::record::data::field::Value<'_>>> {
        let value = Data::get(self, tag)?;
        Some(Ok(value_buf_to_value(value)))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<([u8; 2], crate::alignment::record::data::field::Value<'_>)>>
            + '_,
    > {
        Box::new(Data::iter(self).map(|(tag, value)| {
            let raw_tag = *tag.as_ref();
            Ok((raw_tag, value_buf_to_value(value)))
        }))
    }
}

impl crate::alignment::record::Data for Data {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &[u8; 2],
    ) -> Option<io::Result<crate::alignment::record::data::field::Value<'_>>> {
        let value = self.get(tag)?;
        Some(Ok(value_buf_to_value(value)))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<([u8; 2], crate::alignment::record::data::field::Value<'_>)>>
            + '_,
    > {
        Box::new(self.iter().map(|(tag, value)| {
            let raw_tag = *tag.as_ref();
            Ok((raw_tag, value_buf_to_value(value)))
        }))
    }
}

fn value_buf_to_value(value: &Value) -> crate::alignment::record::data::field::Value<'_> {
    use crate::alignment::record_buf::data::field::value::Array;

    match value {
        Value::Character(c) => crate::alignment::record::data::field::Value::Character(*c),
        Value::Int8(n) => crate::alignment::record::data::field::Value::Int8(*n),
        Value::UInt8(n) => crate::alignment::record::data::field::Value::UInt8(*n),
        Value::Int16(n) => crate::alignment::record::data::field::Value::Int16(*n),
        Value::UInt16(n) => crate::alignment::record::data::field::Value::UInt16(*n),
        Value::Int32(n) => crate::alignment::record::data::field::Value::Int32(*n),
        Value::UInt32(n) => crate::alignment::record::data::field::Value::UInt32(*n),
        Value::Float(n) => crate::alignment::record::data::field::Value::Float(*n),
        Value::String(s) => crate::alignment::record::data::field::Value::String(s),
        Value::Hex(s) => crate::alignment::record::data::field::Value::Hex(s),
        Value::Array(Array::Int8(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::Int8(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::UInt8(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::UInt8(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::Int16(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::Int16(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::UInt16(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::UInt16(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::Int32(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::Int32(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::UInt32(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::UInt32(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
        Value::Array(Array::Float(values)) => crate::alignment::record::data::field::Value::Array(
            crate::alignment::record::data::field::value::Array::Float(Box::new(Values::new(
                values.as_ref(),
            ))),
        ),
    }
}

struct Values<'a, N>(&'a [N]);

impl<'a, N> Values<'a, N>
where
    N: Copy,
{
    fn new(values: &'a [N]) -> Self {
        Self(values)
    }
}

impl<'a, N> crate::alignment::record::data::field::value::array::Values<'a, N> for Values<'a, N>
where
    N: Copy,
{
    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
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
    use super::field::tag;

    use super::*;

    #[test]
    fn test_remove_with_multiple_removes() -> Result<(), field::tag::ParseError> {
        let zz = "zz".parse()?;

        let mut data: Data = [
            (tag::ALIGNMENT_HIT_COUNT, Value::from(2)),
            (tag::EDIT_DISTANCE, Value::from(1)),
            (zz, Value::from(0)),
        ]
        .into_iter()
        .collect();

        data.remove(&tag::EDIT_DISTANCE);
        data.remove(&zz);
        data.remove(&tag::ALIGNMENT_HIT_COUNT);

        assert!(data.is_empty());

        Ok(())
    }

    #[test]
    fn test_from_iterator() {
        let actual: Data = [
            (tag::READ_GROUP, Value::from("rg0")),
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
        ]
        .into_iter()
        .collect();

        let expected: Data = [
            (tag::READ_GROUP, Value::from("rg0")),
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
        ]
        .into_iter()
        .collect();

        assert_eq!(expected, actual);
    }
}
