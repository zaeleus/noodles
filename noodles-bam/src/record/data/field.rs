//! BAM record data field and values.

use std::io;

use noodles_sam::{
    self as sam,
    record::data::field::{Tag, Value},
};

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
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
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
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert_eq!(field.tag(), Tag::AlignmentHitCount);
    /// ```
    pub fn tag(&self) -> Tag {
        self.tag
    }

    /// Returns the data field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert_eq!(field.value(), &Value::Int32(1));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl From<Field> for sam::record::data::Field {
    fn from(field: Field) -> Self {
        Self::new(field.tag, field.value)
    }
}

impl TryFrom<Field> for Vec<u8> {
    type Error = io::Error;

    fn try_from(field: Field) -> Result<Self, Self::Error> {
        use crate::writer::record::data::write_field;

        let mut buf = Vec::new();
        write_field(&mut buf, &field.into())?;
        Ok(buf)
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
            sam::record::data::field::Value::Int32(1),
        );

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_try_from_field_for_vec_u8() -> io::Result<()> {
        let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));

        let actual = <Vec<u8>>::try_from(field)?;
        let expected = vec![b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00];

        assert_eq!(actual, expected);

        Ok(())
    }
}
