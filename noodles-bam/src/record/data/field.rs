//! BAM record data field and values.

pub mod value;

pub use self::value::Value;

use noodles_sam::record::data::field::Tag;

/// A BAM record data field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    tag: Tag,
    value: Value,
}

impl Field {
    /// Creates a data field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::{field::Value, Field};
    /// use noodles_sam::record::data::field::Tag;
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// ```
    pub fn new(tag: Tag, value: Value) -> Self {
        Self { tag, value }
    }

    /// Returns the data field tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::{field::Value, Field};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    ///
    /// assert_eq!(field.tag(), &Tag::AlignmentHitCount);
    /// ```
    pub fn tag(&self) -> &Tag {
        &self.tag
    }

    /// Returns the data field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::{field::Value, Field};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    ///
    /// assert_eq!(field.value(), &Value::Int32(1));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}
