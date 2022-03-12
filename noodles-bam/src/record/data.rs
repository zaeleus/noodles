//! BAM record data and fields.

mod bounds;
pub mod field;
mod fields;

pub(crate) use self::bounds::Bounds;
pub use self::{field::Field, fields::Fields};

use std::{error, fmt, io};

use noodles_sam::{self as sam, record::data::field::Tag};

/// BAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Data {
    data: Vec<u8>,
    bounds: Bounds,
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

    /// Removes all fields from the data map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::Data;
    ///
    /// let mut data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// assert_eq!(data.len(), 2);
    ///
    /// data.clear();
    ///
    /// assert!(data.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn clear(&mut self) {
        self.data.clear();
        self.bounds.clear();
    }

    /// Returns a field by the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::Field, Data};
    /// use noodles_sam::record::data::field::{Tag, Value};
    ///
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// let rg = data.get(Tag::ReadGroup).transpose()?;
    /// assert_eq!(rg, Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0")))));
    ///
    /// assert!(data.get(Tag::AlignmentScore).is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn get(&self, tag: Tag) -> Option<io::Result<Field>> {
        for (i, result) in self.keys().enumerate() {
            match result {
                Ok(t) => {
                    if t == tag {
                        return self.get_index(i);
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    /// Returns a field by an index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::Field, Data};
    /// use noodles_sam::record::data::field::{Tag, Value};
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

    /// Inserts a field into the data map.
    ///
    /// This uses the field tag as the key and field as the value.
    ///
    /// If the tag already exists in the map, the existing field is replaced by the new one, and
    /// the existing field is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::Field, Data};
    /// use noodles_sam::record::data::field::{Tag, Value};
    ///
    /// let mut data = Data::default();
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::UInt8(1));
    /// data.insert(nh.clone()).transpose()?;
    ///
    /// assert_eq!(data.len(), 1);
    /// assert_eq!(data.get_index(0).transpose()?, Some(nh));
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn insert(&mut self, field: Field) -> Option<io::Result<Field>> {
        for (i, result) in self.values().enumerate() {
            match result {
                Ok(f) => {
                    if f.tag() == field.tag() {
                        match self.splice(i, field) {
                            Ok(_) => return Some(Ok(f)),
                            Err(e) => return Some(Err(e)),
                        }
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        if let Err(e) = self.push(field) {
            Some(Err(e))
        } else {
            None
        }
    }

    fn splice(&mut self, i: usize, field: Field) -> io::Result<()> {
        let range = self.bounds.get(i).expect("index out of bounds");
        let buf = <Vec<u8>>::try_from(field)?;
        self.data.splice(range, buf);

        self.index()?;

        Ok(())
    }

    fn push(&mut self, field: Field) -> io::Result<()> {
        let buf = <Vec<u8>>::try_from(field)?;
        self.data.extend_from_slice(&buf);
        self.bounds.push(self.data.len());
        Ok(())
    }

    /// Returns an iterator over data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::Field, Data};
    /// use noodles_sam::record::data::field::{Tag, Value};
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

    /// Returns an iterator over all tags.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::Data;
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    ///     b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
    /// ])?;
    ///
    /// let mut keys = data.keys();
    ///
    /// assert_eq!(keys.next().transpose()?, Some(Tag::AlignmentHitCount));
    /// assert_eq!(keys.next().transpose()?, Some(Tag::ReadGroup));
    /// assert!(keys.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = io::Result<Tag>> + '_ {
        let mut start = 0;

        self.bounds.as_ref().iter().map(move |&end| {
            let buf = &self.data[start..end];
            let raw_tag = [buf[0], buf[1]];

            start = end;

            Tag::try_from(raw_tag).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    /// Returns an iterator over all fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::Field, Data};
    /// use noodles_sam::record::data::field::{Tag, Value};
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
        self.bounds.update(&self.data[..])
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

impl TryFrom<Vec<Field>> for Data {
    type Error = io::Error;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        let mut buf = Vec::new();

        for field in fields {
            let raw_field = <Vec<u8>>::try_from(field)?;
            buf.extend_from_slice(&raw_field);
        }

        Self::try_from(buf)
    }
}

/// An error returned when BAM data fails to convert to SAM data.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromDataError {
    /// A field is invalid.
    InvalidField,
    /// The data is invalid.
    DuplicateTag(Tag),
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
    use sam::record::data::field::Value;

    use super::*;

    #[test]
    fn test_insert() -> io::Result<()> {
        let mut data = Data::try_from(vec![
            b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ])?;

        // push
        let as_ = Field::new(Tag::AlignmentScore, Value::UInt8(13));
        assert!(data.insert(as_).transpose()?.is_none());
        assert_eq!(
            data.as_ref(),
            [
                b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
                b'A', b'S', b'C', 0x0d, // AS:C:13
            ]
        );

        // replace
        let nh = Field::new(Tag::AlignmentHitCount, Value::UInt8(2));
        let expected_nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
        assert_eq!(data.insert(nh).transpose()?, Some(expected_nh));
        assert_eq!(
            data.as_ref(),
            [
                b'N', b'H', b'C', 0x02, // NH:C:2
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
                b'A', b'S', b'C', 0x0d, // AS:C:13
            ]
        );

        Ok(())
    }

    #[test]
    fn test_try_from_vec_field_for_data() -> io::Result<()> {
        let fields = vec![
            Field::new(Tag::AlignmentHitCount, Value::UInt8(1)),
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ];

        let data = Data::try_from(fields)?;

        assert_eq!(
            data.as_ref(),
            vec![
                b'N', b'H', b'C', 0x01, // NH:i:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
        );

        Ok(())
    }

    #[test]
    fn test_try_from_data_for_sam_record_data() -> Result<(), Box<dyn std::error::Error>> {
        let data = Data::try_from(vec![
            b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ])?;

        let actual = sam::record::Data::try_from(&data)?;
        let expected = sam::record::Data::try_from(vec![
            sam::record::data::Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
            sam::record::data::Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
        ])?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
