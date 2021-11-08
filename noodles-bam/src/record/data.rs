//! BAM record data and fields.

mod bounds;
pub mod field;
mod fields;

pub(crate) use self::bounds::Bounds;
pub use self::{field::Field, fields::Fields};

use std::{error, fmt, io};

use noodles_sam as sam;

/// BAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Data {
    data: Vec<u8>,
    pub(crate) bounds: Bounds,
}

impl Data {
    /// Returns the number of data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Data;
    /// let data = Data::default();
    /// assert_eq!(data.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.bounds.len()
    }

    /// Returns whether there are any data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Data;
    /// let data = Data::default();
    /// assert!(data.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns a field by an index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::{field::Value, Field}, Data};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// let rg = data.get_index(1).transpose()?;
    /// assert_eq!(rg, Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0")))));
    ///
    /// assert!(data.get_index(2).is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn get_index(&self, i: usize) -> Option<io::Result<Field>> {
        use crate::reader::record::data::read_field;

        self.bounds.get(i).map(|range| {
            let mut reader = &self.data[range];
            read_field(&mut reader).transpose().unwrap()
        })
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
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// let mut fields = data.values();
    ///
    /// assert_eq!(
    ///     fields.next().transpose()?,
    ///     Some(Field::new(Tag::AlignmentHitCount, Value::Int32(1)))
    /// );
    ///
    /// assert_eq!(
    ///     fields.next().transpose()?,
    ///     Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))))
    /// );
    ///
    /// assert!(fields.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    #[deprecated(since = "0.8.0", note = "Use `Data::values` instead.")]
    pub fn fields(&self) -> Fields<&[u8]> {
        self.values()
    }

    /// Returns an iterator over all fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::{field::Value, Field}, Data};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// let mut values = data.values();
    ///
    /// assert_eq!(
    ///     values.next().transpose()?,
    ///     Some(Field::new(Tag::AlignmentHitCount, Value::Int32(1)))
    /// );
    ///
    /// assert_eq!(
    ///     values.next().transpose()?,
    ///     Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))))
    /// );
    ///
    /// assert!(values.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn values(&self) -> Fields<&[u8]> {
        Fields::new(&self.data)
    }

    pub(crate) fn index(&mut self) -> io::Result<()> {
        self.bounds = Bounds::try_from_buf(&self.data[..])?;
        Ok(())
    }
}

impl AsRef<[u8]> for Data {
    fn as_ref(&self) -> &[u8] {
        &self.data
    }
}

impl AsMut<Vec<u8>> for Data {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }
}

impl fmt::Debug for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.values()).finish()
    }
}

impl TryFrom<Vec<u8>> for Data {
    type Error = io::Error;

    fn try_from(data: Vec<u8>) -> Result<Self, Self::Error> {
        let mut data = Self {
            data,
            bounds: Bounds::default(),
        };

        data.index()?;

        Ok(data)
    }
}

/// An error returned when BAM data fails to convert to SAM data.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromDataError {
    /// A field is invalid.
    InvalidField,
    /// The data is invalid.
    DuplicateTag(sam::record::data::field::Tag),
}

impl error::Error for TryFromDataError {}

impl fmt::Display for TryFromDataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField => write!(f, "invalid field"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
    }
}

impl TryFrom<&Data> for sam::record::Data {
    type Error = TryFromDataError;

    fn try_from(data: &Data) -> Result<Self, Self::Error> {
        let mut sam_data = Self::default();

        for result in data.values() {
            let field = result.map_err(|_| TryFromDataError::InvalidField)?;
            let sam_field = sam::record::data::Field::from(field);

            if let Some(f) = sam_data.insert(sam_field) {
                return Err(TryFromDataError::DuplicateTag(f.tag()));
            }
        }

        Ok(sam_data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_data_for_sam_record_data() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::data::{
            field::{Tag, Value},
            Field,
        };

        let data = Data::try_from(vec![
            b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ])?;

        let actual = sam::record::Data::try_from(&data)?;
        let expected = sam::record::Data::try_from(vec![
            Field::new(Tag::AlignmentHitCount, Value::Int(1)),
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ])?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
