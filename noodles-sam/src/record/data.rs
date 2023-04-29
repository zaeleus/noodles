//! SAM record data and fields.

pub mod field;

use std::{
    error,
    fmt::{self, Write},
    mem,
    str::FromStr,
};

const DELIMITER: char = '\t';

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Data {
    fields: Vec<(field::Tag, field::Value)>,
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
    pub fn get<K>(&self, tag: &K) -> Option<&field::Value>
    where
        K: indexmap::Equivalent<field::Tag>,
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
        K: indexmap::Equivalent<field::Tag>,
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
    pub fn iter(&self) -> impl Iterator<Item = (field::Tag, &field::Value)> {
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
    pub fn keys(&self) -> impl Iterator<Item = field::Tag> + '_ {
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
    pub fn values(&self) -> impl Iterator<Item = &field::Value> {
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
    pub fn insert(
        &mut self,
        tag: field::Tag,
        value: field::Value,
    ) -> Option<(field::Tag, field::Value)> {
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
    /// let rg = (tag::READ_GROUP, Value::String(String::from("rg0")));
    /// let md = (tag::ALIGNMENT_SCORE, Value::from(98));
    /// let mut data: Data = [nh.clone(), rg.clone(), md.clone()].into_iter().collect();
    ///
    /// assert_eq!(data.remove(&tag::ALIGNMENT_HIT_COUNT), Some(nh));
    /// assert!(data.remove(&tag::COMMENT).is_none());
    ///
    /// let expected = [md, rg].into_iter().collect();
    /// assert_eq!(data, expected);
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn remove<K>(&mut self, tag: &K) -> Option<(field::Tag, field::Value)>
    where
        K: indexmap::Equivalent<field::Tag>,
    {
        self.swap_remove(tag)
    }

    fn swap_remove<K>(&mut self, tag: &K) -> Option<(field::Tag, field::Value)>
    where
        K: indexmap::Equivalent<field::Tag>,
    {
        self.get_index_of(tag).map(|i| self.fields.swap_remove(i))
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use field::Type;

        for (i, (tag, value)) in self.iter().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            let ty = if value.is_int() {
                Type::Int32
            } else {
                value.ty()
            };

            write!(f, "{tag}:{ty}:{value}")?;
        }

        Ok(())
    }
}

impl Extend<(field::Tag, field::Value)> for Data {
    fn extend<T: IntoIterator<Item = (field::Tag, field::Value)>>(&mut self, iter: T) {
        for (tag, value) in iter {
            self.insert(tag, value);
        }
    }
}

/// An error returned when raw SAM record data fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input data contains an invalid field.
    InvalidField(field::ParseError),
    /// A tag is duplicated.
    ///
    /// § 1.5 The alignment section: optional fields (2021-01-07): "Each `TAG` can only appear once
    /// in one alignment line."
    DuplicateTag(field::Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::DuplicateTag(_) => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => f.write_str("invalid field"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

impl FromIterator<(field::Tag, field::Value)> for Data {
    fn from_iter<T: IntoIterator<Item = (field::Tag, field::Value)>>(iter: T) -> Self {
        let mut data = Self::default();
        data.extend(iter);
        data
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        let mut data = Self::default();

        for s in s.split(DELIMITER) {
            let (tag, value) = parse_field(s).map_err(ParseError::InvalidField)?;

            if data.insert(tag, value).is_some() {
                return Err(ParseError::DuplicateTag(tag));
            }
        }

        Ok(data)
    }
}

fn parse_field(s: &str) -> Result<(field::Tag, field::Value), field::ParseError> {
    let (raw_tag, rest) = s.split_once(':').ok_or(field::ParseError::Invalid)?;
    let tag = raw_tag.parse().map_err(field::ParseError::InvalidTag)?;

    let (raw_ty, raw_value) = rest.split_once(':').ok_or(field::ParseError::Invalid)?;
    let ty = raw_ty.parse().map_err(field::ParseError::InvalidType)?;
    let value =
        field::Value::from_str_type(raw_value, ty).map_err(|_| field::ParseError::InvalidValue)?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::field::{tag, Value};

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
    fn test_fmt() {
        let data: Data = [
            (tag::READ_GROUP, Value::String(String::from("rg0"))),
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
        ]
        .into_iter()
        .collect();

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);
    }

    #[test]
    fn test_from_iterator() {
        let actual: Data = [
            (tag::READ_GROUP, Value::String(String::from("rg0"))),
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
        ]
        .into_iter()
        .collect();

        let expected: Data = [
            (tag::READ_GROUP, Value::String(String::from("rg0"))),
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
        ]
        .into_iter()
        .collect();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_from_str() {
        assert_eq!("".parse(), Ok(Data::default()));

        assert_eq!(
            "RG:Z:rg0\tNH:i:1".parse::<Data>(),
            Ok([
                (tag::READ_GROUP, Value::String(String::from("rg0"))),
                (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            ]
            .into_iter()
            .collect())
        );

        assert_eq!(
            "NH:i:1\tNH:i:1".parse::<Data>(),
            Err(ParseError::DuplicateTag(tag::ALIGNMENT_HIT_COUNT))
        );
    }
}
