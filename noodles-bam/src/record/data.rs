//! BAM record data and fields.

pub mod field;
pub mod reader;

pub use self::{field::Field, reader::Reader};

use std::{convert::TryFrom, io, ops::Deref};

use noodles_sam as sam;

use self::reader::Fields;

/// BAM record data.
///
/// This is also called optional fields.
#[derive(Debug)]
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    /// Creates data from raw data data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Data;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let raw_data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let data = Data::new(&raw_data);
    /// ```
    pub fn new(bytes: &[u8]) -> Data<'_> {
        Data(bytes)
    }

    /// Returns an iterator over data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::{field::Value, Field}, Data};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let raw_data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let data = Data::new(&raw_data);
    ///
    /// let mut fields = data.fields();
    ///
    /// let field = fields.next().unwrap()?;
    /// assert_eq!(field, Field::new(Tag::AlignmentHitCount, Value::Int32(1)));
    ///
    /// let field = fields.next().unwrap()?;
    /// assert_eq!(field, Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))));
    ///
    /// assert!(fields.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn fields(&self) -> Fields<&[u8]> {
        let reader = Reader::new(self.0);
        reader.fields()
    }
}

impl<'a> Deref for Data<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.0
    }
}

impl<'a> TryFrom<Data<'a>> for sam::record::Data {
    type Error = io::Error;

    fn try_from(data: Data<'_>) -> Result<Self, Self::Error> {
        use field::Value as BamValue;
        use sam::record::data::field::Value as SamValue;

        let mut sam_fields = Vec::new();

        for result in data.fields() {
            let field = result?;
            let tag = field.tag();

            let value = match field.value() {
                BamValue::Char(c) => SamValue::Char(*c),
                BamValue::Int8(n) => SamValue::Int32(*n as i32),
                BamValue::UInt8(n) => SamValue::Int32(*n as i32),
                BamValue::Int16(n) => SamValue::Int32(*n as i32),
                BamValue::UInt16(n) => SamValue::Int32(*n as i32),
                BamValue::Int32(n) => SamValue::Int32(*n),
                // FIXME: lossy conversion
                BamValue::UInt32(n) => SamValue::Int32(*n as i32),
                BamValue::Float(n) => SamValue::Float(*n),
                BamValue::String(s) => SamValue::String(s.clone()),
                BamValue::Hex(s) => SamValue::Hex(s.clone()),
                BamValue::Int8Array(values) => SamValue::Int8Array(values.clone()),
                BamValue::UInt8Array(values) => SamValue::UInt8Array(values.clone()),
                BamValue::Int16Array(values) => SamValue::Int16Array(values.clone()),
                BamValue::UInt16Array(values) => SamValue::UInt16Array(values.clone()),
                BamValue::Int32Array(values) => SamValue::Int32Array(values.clone()),
                BamValue::UInt32Array(values) => SamValue::UInt32Array(values.clone()),
                BamValue::FloatArray(values) => SamValue::FloatArray(values.clone()),
            };

            let sam_field = sam::record::data::Field::new(tag.clone(), value);
            sam_fields.push(sam_field);
        }

        Ok(sam::record::Data::from(sam_fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_data_for_sam_record_data() -> io::Result<()> {
        use sam::record::data::{
            field::{Tag, Value},
            Field,
        };

        let raw_data = [
            0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00, // NH:i:1
            0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00, // RG:Z:rg0
        ];
        let data = Data::new(&raw_data);

        let actual = sam::record::Data::try_from(data)?;
        let expected = sam::record::Data::from(vec![
            Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }
}
