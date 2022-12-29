//! SAM record data and fields.

pub mod field;

pub use self::field::Field;

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
#[derive(Default, Clone, PartialEq)]
pub struct Data {
    fields: Vec<Field>,
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
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let mut data = Data::try_from(vec![
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    /// ])?;
    ///
    /// assert_eq!(data.len(), 1);
    ///
    /// data.clear();
    ///
    /// assert!(data.is_empty());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn clear(&mut self) {
        self.fields.clear();
    }

    /// Returns a reference to the value of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh.clone()])?;
    ///
    /// assert_eq!(data.get(Tag::AlignmentHitCount), Some(nh.value()));
    /// assert!(data.get(Tag::ReadGroup).is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn get(&self, tag: field::Tag) -> Option<&field::Value> {
        self.fields
            .iter()
            .find(|f| f.tag() == tag)
            .map(|f| f.value())
    }

    /// Returns the index of the field of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh])?;
    ///
    /// assert_eq!(data.get_index_of(Tag::AlignmentHitCount), Some(0));
    /// assert!(data.get_index_of(Tag::ReadGroup).is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn get_index_of(&self, tag: field::Tag) -> Option<usize> {
        self.fields.iter().position(|f| f.tag() == tag)
    }

    /// Returns an iterator over all tag-value pairs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh.clone()])?;
    ///
    /// let mut fields = data.iter();
    /// assert_eq!(fields.next(), Some((nh.tag(), nh.value())));
    /// assert!(fields.next().is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (field::Tag, &field::Value)> {
        self.fields.iter().map(|f| (f.tag(), f.value()))
    }

    /// Returns an iterator over all tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh])?;
    ///
    /// let mut keys = data.keys();
    /// assert_eq!(keys.next(), Some(Tag::AlignmentHitCount));
    /// assert!(keys.next().is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = field::Tag> + '_ {
        self.fields.iter().map(|field| field.tag())
    }

    /// Returns an iterator over all values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh.clone()])?;
    ///
    /// let mut values = data.values();
    /// assert_eq!(values.next(), Some(nh.value()));
    /// assert!(values.next().is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &field::Value> {
        self.fields.iter().map(|field| field.value())
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
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    /// let mut data = Data::default();
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// data.insert(nh);
    /// ```
    pub fn insert(&mut self, field: Field) -> Option<Field> {
        match self.get_index_of(field.tag()) {
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
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let rg = Field::new(Tag::ReadGroup, Value::String(String::from("rg0")));
    /// let md = Field::new(Tag::AlignmentScore, Value::Int32(98));
    /// let mut data = Data::try_from(vec![nh.clone(), rg.clone(), md.clone()])?;
    ///
    /// assert_eq!(data.remove(Tag::AlignmentHitCount), Some(nh));
    /// assert!(data.remove(Tag::Comment).is_none());
    ///
    /// let expected = Data::try_from(vec![md, rg])?;
    /// assert_eq!(data, expected);
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn remove(&mut self, tag: field::Tag) -> Option<Field> {
        self.swap_remove(tag)
    }

    fn swap_remove(&mut self, tag: field::Tag) -> Option<Field> {
        self.get_index_of(tag).map(|i| self.fields.swap_remove(i))
    }
}

impl fmt::Debug for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.fields.iter()).finish()
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use field::value::Type;

        for (i, (tag, value)) in self.iter().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            let ty = if value.is_int() {
                Type::Int32
            } else {
                value.ty()
            };

            write!(f, "{}:{}:{}", tag, ty, value)?;
        }

        Ok(())
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
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
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
            let field = s.parse().map_err(ParseError::InvalidField)?;

            if let Some(f) = data.insert(field) {
                return Err(ParseError::DuplicateTag(f.tag()));
            }
        }

        Ok(data)
    }
}

impl TryFrom<Vec<Field>> for Data {
    type Error = ParseError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        let mut data = Self::default();

        for field in fields {
            if let Some(f) = data.insert(field) {
                return Err(ParseError::DuplicateTag(f.tag()));
            }
        }

        Ok(data)
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_remove_with_multiple_removes() -> Result<(), Box<dyn std::error::Error>> {
        let zz = "zz".parse()?;

        let mut data = Data::try_from(vec![
            Field::new(Tag::AlignmentHitCount, Value::UInt8(2)),
            Field::new(Tag::EditDistance, Value::UInt8(1)),
            Field::new(zz, Value::UInt8(0)),
        ])?;

        data.remove(Tag::EditDistance);
        data.remove(zz);
        data.remove(Tag::AlignmentHitCount);

        assert!(data.is_empty());

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let data = Data::try_from(vec![
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
            Field::new(Tag::AlignmentHitCount, Value::UInt8(1)),
        ])?;

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("".parse(), Ok(Data::default()));

        assert_eq!(
            "RG:Z:rg0\tNH:i:1".parse(),
            Ok(Data::try_from(vec![
                Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
                Field::new(Tag::AlignmentHitCount, Value::UInt8(1)),
            ])?)
        );

        assert_eq!(
            "NH:i:1\tNH:i:1".parse::<Data>(),
            Err(ParseError::DuplicateTag(Tag::AlignmentHitCount))
        );

        Ok(())
    }
}
