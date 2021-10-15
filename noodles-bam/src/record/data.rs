//! BAM record data and fields.

pub mod field;
mod fields;

pub use self::{field::Field, fields::Fields};

use std::{convert::TryFrom, error, fmt, ops::Deref};

use noodles_sam as sam;

/// BAM record data.
///
/// This is also called optional fields.
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    /// Creates data from raw data data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Data;
    ///
    /// let raw_data = [
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
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
    /// let raw_data = [
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ];
    /// let data = Data::new(&raw_data);
    ///
    /// let mut fields = data.fields();
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
    pub fn fields(&self) -> Fields<&[u8]> {
        Fields::new(self.0)
    }
}

impl<'a> Deref for Data<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.0
    }
}

impl<'a> fmt::Debug for Data<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.fields()).finish()
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

impl<'a> TryFrom<Data<'a>> for sam::record::Data {
    type Error = TryFromDataError;

    fn try_from(data: Data<'_>) -> Result<Self, Self::Error> {
        let mut sam_data = Self::default();

        for result in data.fields() {
            let field = result.map_err(|_| TryFromDataError::InvalidField)?;

            let tag = field.tag();
            let sam_field = sam::record::data::Field::from(field);

            if sam_data.insert(tag, sam_field).is_some() {
                return Err(TryFromDataError::DuplicateTag(tag));
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

        let raw_data = [
            b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ];
        let data = Data::new(&raw_data);

        let actual = sam::record::Data::try_from(data)?;
        let expected = sam::record::Data::try_from(vec![
            Field::new(Tag::AlignmentHitCount, Value::Int(1)),
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ])?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
