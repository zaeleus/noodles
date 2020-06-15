//! BAM record data and fields.

pub mod field;
pub mod reader;

pub use self::{field::Field, reader::Reader};

use std::ops::Deref;

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
