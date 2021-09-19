//! BAM record data reader and iterators.

mod fields;

pub use self::fields::Fields;

use std::io::{self, BufRead};

use super::{
    field::{value::Type, Value},
    Field,
};

/// A BAM record data reader.
pub struct Reader<R>
where
    R: BufRead,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a BAM record data reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Reader;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let reader = Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a BAM record data field value.
    ///
    /// The stream is expected to be at the start of the value, i.e., after the tag and data type.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::data::{field::{value::Type, Value}, Reader};
    ///
    /// let data = [0x01, 0x00, 0x00, 0x00];
    /// let mut reader = Reader::new(&data[..]);
    ///
    /// let value = reader.read_value_type(Type::Int32)?;
    ///
    /// assert_eq!(value, Value::Int32(1));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_value_type(&mut self, ty: Type) -> io::Result<Value> {
        use crate::reader::record::data::field::read_value;

        read_value(&mut self.inner, ty)
    }

    /// Returns an iterator over data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::{field::Value, Field, Reader}};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let reader = Reader::new(&data[..]);
    ///
    /// let mut fields = reader.fields();
    ///
    /// let field = fields.next().transpose()?;
    /// assert_eq!(field, Some(Field::new(Tag::AlignmentHitCount, Value::Int32(1))));
    ///
    /// let field = fields.next().transpose()?;
    /// assert_eq!(field, Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0")))));
    ///
    /// assert!(fields.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn fields(self) -> Fields<R> {
        Fields::new(self)
    }

    fn read_field(&mut self) -> io::Result<Option<Field>> {
        use crate::reader::record::data::read_field;

        match read_field(&mut self.inner) {
            Ok(field) => Ok(Some(field)),
            Err(_) => Ok(None),
        }
    }
}
