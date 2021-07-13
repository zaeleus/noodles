//! CRAM record tag and fields.

mod key;

pub use self::key::Key;

use noodles_bam::record::data::field::Value;

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
    /// use noodles_bam::record::data::field::Value;
    /// use noodles_cram::record::{tag::Key, Tag};
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
    /// use noodles_bam::record::data::field::Value;
    /// use noodles_cram::record::{tag::Key, Tag};
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
    /// use noodles_bam::record::data::field::Value;
    /// use noodles_cram::record::{tag::Key, Tag};
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
