//! CRAM record tag and fields.

use std::{error, fmt};

mod key;

pub use self::key::Key;

use noodles_sam::{self as sam, record::data::field::Value};

/// A CRAM record tag.
#[derive(Clone, Debug, PartialEq)]
pub struct Tag {
    key: Key,
    value: Value,
}

impl Tag {
    /// Creates a CRAM record tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::{tag::Key, Tag};
    /// use noodles_sam::record::data::field::Value;
    ///
    /// let value = Value::Int8(1);
    /// let key = Key::new([b'N', b'H'], value.ty());
    ///
    /// let tag = Tag::new(key, value);
    /// ```
    pub fn new(key: Key, value: Value) -> Self {
        Self { key, value }
    }

    /// Returns the tag key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::{tag::Key, Tag};
    /// use noodles_sam::record::data::field::Value;
    ///
    /// let value = Value::Int8(1);
    /// let key = Key::new([b'N', b'H'], value.ty());
    ///
    /// let tag = Tag::new(key.clone(), value);
    /// assert_eq!(tag.key(), key);
    /// ```
    pub fn key(&self) -> Key {
        self.key
    }

    /// Returns the tag value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::{tag::Key, Tag};
    /// use noodles_sam::record::data::field::Value;
    ///
    /// let value = Value::Int8(1);
    /// let key = Key::new([b'N', b'H'], value.ty());
    ///
    /// let tag = Tag::new(key, value.clone());
    /// assert_eq!(tag.value(), &value);
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

/// An error returned when a CRAM record tag fails to convert to a SAM record data field.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromTagError {
    /// The tag is invalid.
    InvalidTag,
}

impl error::Error for TryFromTagError {}

impl fmt::Display for TryFromTagError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidTag => f.write_str("invalid tag"),
        }
    }
}

impl TryFrom<Tag> for sam::record::data::Field {
    type Error = TryFromTagError;

    fn try_from(tag: Tag) -> Result<Self, Self::Error> {
        let sam_tag = tag
            .key()
            .tag()
            .try_into()
            .map_err(|_| TryFromTagError::InvalidTag)?;

        Ok(Self::new(sam_tag, tag.value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_tag_for_sam_record_data_field() {
        use sam::record::data::Field;

        let value = Value::Int8(1);
        let key = Key::new([b'N', b'H'], value.ty());
        let tag = Tag::new(key, value);
        let actual = Field::try_from(tag);

        let expected = Ok(Field::new(
            sam::record::data::field::Tag::AlignmentHitCount,
            sam::record::data::field::Value::Int8(1),
        ));

        assert_eq!(actual, expected);
    }
}
