//! BAM record data field and values.

pub mod value;

pub use self::value::Value;

use noodles_sam::{self as sam, record::data::field::Tag};

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

impl From<Field> for sam::record::data::Field {
    fn from(field: Field) -> Self {
        Self::new(field.tag, field.value.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_field_for_sam_record_data_field() {
        let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));

        let actual = sam::record::data::Field::from(field);
        let expected = sam::record::data::Field::new(
            Tag::AlignmentHitCount,
            sam::record::data::field::Value::Int(1),
        );

        assert_eq!(actual, expected);
    }
}
