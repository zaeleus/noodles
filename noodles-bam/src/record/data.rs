//! BAM record data and fields.

pub mod field;
mod fields;

pub use self::{field::Field, fields::Fields};

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
};

use noodles_sam as sam;

/// BAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Data(Vec<u8>);

impl Data {
    /// Creates data from raw data data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Data;
    ///
    /// let data = Data::from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ]);
    /// ```
    #[deprecated(since = "0.8.0", note = "Use `Data::from::<Vec<u8>>` instead.")]
    pub fn new(data: Vec<u8>) -> Data {
        Data::from(data)
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
    /// let data = Data::from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ]);
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
        Fields::new(self)
    }
}

impl Deref for Data {
    type Target = Vec<u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Data {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Debug for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.fields()).finish()
    }
}

impl From<Vec<u8>> for Data {
    fn from(data: Vec<u8>) -> Self {
        Self(data)
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

        for result in data.fields() {
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

        let data = Data::from(vec![
            b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ]);

        let actual = sam::record::Data::try_from(&data)?;
        let expected = sam::record::Data::try_from(vec![
            Field::new(Tag::AlignmentHitCount, Value::Int(1)),
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ])?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
