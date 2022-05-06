//! SAM record data and fields.

pub mod field;

pub use self::field::Field;

use std::{
    error,
    fmt::{self, Write},
    mem,
    num::NonZeroU16,
    str::FromStr,
};

use rustc_hash::FxHashMap;

type StandardFieldIndices = [Option<NonZeroU16>; 55];
type OtherFieldIndices = FxHashMap<field::Tag, u16>;

const DELIMITER: char = '\t';

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Clone, PartialEq)]
pub struct Data {
    standard_field_indices: StandardFieldIndices,
    other_field_indices: OtherFieldIndices,
    fields: Vec<Field>,
}

impl Data {
    /// Returns the number of data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Data;
    /// let data = Data::default();
    /// assert_eq!(data.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// Returns whether there are any data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Data;
    /// let data = Data::default();
    /// assert!(data.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Removes all data fields from the data map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let mut data = Data::try_from(vec![
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    /// ])?;
    ///
    /// assert_eq!(data.len(), 1);
    ///
    /// data.clear();
    ///
    /// assert!(data.is_empty());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn clear(&mut self) {
        self.standard_field_indices.fill(None);
        self.other_field_indices.clear();
        self.fields.clear();
    }

    /// Returns a reference to the field of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh.clone()])?;
    ///
    /// assert_eq!(data.get(Tag::AlignmentHitCount), Some(&nh));
    /// assert!(data.get(Tag::ReadGroup).is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn get(&self, tag: field::Tag) -> Option<&Field> {
        self.get_index_of(tag).and_then(|j| self.fields.get(j))
    }

    /// Returns the index of the field of the given tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh])?;
    ///
    /// assert_eq!(data.get_index_of(Tag::AlignmentHitCount), Some(0));
    /// assert!(data.get_index_of(Tag::ReadGroup).is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn get_index_of(&self, tag: field::Tag) -> Option<usize> {
        match tag_to_index(tag) {
            Some(i) => self.get_normalized_standard_field_index(i),
            None => self.get_normalized_other_field_index(tag),
        }
    }

    /// Returns an iterator over all tags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh])?;
    ///
    /// let mut keys = data.keys();
    /// assert_eq!(keys.next(), Some(Tag::AlignmentHitCount));
    /// assert!(keys.next().is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = field::Tag> + '_ {
        self.fields.iter().map(|field| field.tag())
    }

    /// Returns an iterator over all fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let data = Data::try_from(vec![nh.clone()])?;
    ///
    /// let mut values = data.values();
    /// assert_eq!(values.next(), Some(&nh));
    /// assert!(values.next().is_none());
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &Field> {
        self.fields.iter()
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
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    /// let mut data = Data::default();
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// data.insert(nh);
    /// ```
    pub fn insert(&mut self, field: Field) -> Option<Field> {
        match self.get_index_of(field.tag()) {
            Some(i) => Some(mem::replace(&mut self.fields[i], field)),
            None => {
                self.push(field);
                None
            }
        }
    }

    /// Removes the field with the given tag.
    ///
    /// The field is returned if it exists.
    ///
    /// This works like [`Vec::swap_remove`]; it does not preserve the order but has a constant
    /// time complexity.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// let rg = Field::new(Tag::ReadGroup, Value::String(String::from("rg0")));
    /// let md = Field::new(Tag::AlignmentScore, Value::Int32(98));
    /// let mut data = Data::try_from(vec![nh.clone(), rg.clone(), md.clone()])?;
    ///
    /// assert_eq!(data.remove(Tag::AlignmentHitCount), Some(nh));
    /// assert!(data.remove(Tag::Comment).is_none());
    ///
    /// let expected = Data::try_from(vec![md, rg])?;
    /// assert_eq!(data, expected);
    /// # Ok::<_, noodles_sam::record::data::ParseError>(())
    /// ```
    pub fn remove(&mut self, tag: field::Tag) -> Option<Field> {
        self.swap_remove(tag)
    }

    fn get_normalized_standard_field_index(&self, i: usize) -> Option<usize> {
        self.standard_field_indices[i].map(|j| {
            // SAFETY: `j` is guaranteed > 0.
            usize::from(u16::from(j) - 1)
        })
    }

    fn get_normalized_other_field_index(&self, tag: field::Tag) -> Option<usize> {
        self.other_field_indices.get(&tag).copied().map(usize::from)
    }

    fn push(&mut self, field: Field) {
        set_index(
            &mut self.standard_field_indices,
            &mut self.other_field_indices,
            field.tag(),
            self.fields.len(),
        );

        self.fields.push(field);
    }

    fn swap_remove(&mut self, tag: field::Tag) -> Option<Field> {
        let i = match tag_to_index(tag) {
            Some(i) => {
                let j = self.get_normalized_standard_field_index(i);
                self.standard_field_indices[i].take();
                j
            }
            None => {
                let j = self.get_normalized_other_field_index(tag);
                self.other_field_indices.remove(&tag);
                j
            }
        };

        i.map(|j| {
            let removed_field = self.fields.swap_remove(j);

            let swapped_field = &self.fields[j];
            set_index(
                &mut self.standard_field_indices,
                &mut self.other_field_indices,
                swapped_field.tag(),
                j,
            );

            removed_field
        })
    }
}

impl Default for Data {
    fn default() -> Self {
        let standard_field_indices = [
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
            None, None, None, None, None, None, None, None, None, None, None, None, None,
        ];

        Self {
            standard_field_indices,
            other_field_indices: FxHashMap::default(),
            fields: Vec::with_capacity(16),
        }
    }
}

impl fmt::Debug for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.fields.iter()).finish()
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.values().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

/// An error returned when raw SAM record data fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input data contains an invalid field.
    InvalidField(field::ParseError),
    /// A tag is duplicated.
    ///
    /// § 1.5 The alignment section: optional fields (2021-01-07): "Each `TAG` can only appear once
    /// in one alignment line."
    DuplicateTag(field::Tag),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => write!(f, "invalid field: {}", e),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        let mut data = Self::default();

        for s in s.split(DELIMITER) {
            let field = s.parse().map_err(ParseError::InvalidField)?;

            if let Some(f) = data.insert(field) {
                return Err(ParseError::DuplicateTag(f.tag()));
            }
        }

        Ok(data)
    }
}

impl TryFrom<Vec<Field>> for Data {
    type Error = ParseError;

    fn try_from(fields: Vec<Field>) -> Result<Self, Self::Error> {
        let mut data = Self::default();

        for field in fields {
            if let Some(f) = data.insert(field) {
                return Err(ParseError::DuplicateTag(f.tag()));
            }
        }

        Ok(data)
    }
}

fn set_index(
    standard_field_indices: &mut StandardFieldIndices,
    other_field_indices: &mut OtherFieldIndices,
    tag: field::Tag,
    i: usize,
) {
    // SAFETY: `i` is guaranteed to be < 3224.
    let j = i as u16;

    match tag_to_index(tag) {
        Some(i) => {
            let j = NonZeroU16::new(j + 1);
            standard_field_indices[i] = j;
        }
        None => {
            other_field_indices.insert(tag, j);
        }
    }
}

fn tag_to_index(tag: field::Tag) -> Option<usize> {
    use field::Tag;

    match tag {
        Tag::MinMappingQuality => Some(0),
        Tag::AlignmentScore => Some(1),
        Tag::SampleBarcodeSequence => Some(2),
        Tag::BaseAlignmentQualityOffsets => Some(3),
        Tag::OriginalUmiQualityScores => Some(4),
        Tag::CellBarcodeId => Some(5),
        Tag::NextHitReferenceSequenceName => Some(6),
        Tag::Cigar => Some(7),
        Tag::ColorEditDistance => Some(8),
        Tag::Comment => Some(9),
        Tag::NextHitPosition => Some(10),
        Tag::ColarQualityScores => Some(11),
        Tag::CellBarcodeSequence => Some(12),
        Tag::ColorSequence => Some(13),
        Tag::CompleteReadAnnotations => Some(14),
        Tag::CellBarcodeQualityScores => Some(15),
        Tag::NextHitSequence => Some(16),
        Tag::SegmentIndex => Some(17),
        Tag::SegmentSuffix => Some(18),
        Tag::AlternativeSequence => Some(19),
        Tag::PerfectHitCount => Some(20),
        Tag::OneDifferenceHitCount => Some(21),
        Tag::TwoDifferenceHitCount => Some(22),
        Tag::HitIndex => Some(23),
        Tag::TotalHitCount => Some(24),
        Tag::Library => Some(25),
        Tag::MateCigar => Some(26),
        Tag::MismatchedPositions => Some(27),
        Tag::UmiId => Some(28),
        Tag::BaseModificationProbabilities => Some(29),
        Tag::BaseModifications => Some(30),
        Tag::MateMappingQuality => Some(31),
        Tag::AlignmentHitCount => Some(32),
        Tag::EditDistance => Some(33),
        Tag::OriginalAlignment => Some(34),
        Tag::OriginalCigar => Some(35),
        Tag::OriginalPosition => Some(36),
        Tag::OriginalQualityScores => Some(37),
        Tag::OriginalUmiBarcodeSequence => Some(38),
        Tag::Program => Some(39),
        Tag::TemplateLikelihood => Some(40),
        Tag::PaddedReadAnnotations => Some(41),
        Tag::PlatformUnit => Some(42),
        Tag::MateQualityScores => Some(43),
        Tag::SampleBarcodeQualityScores => Some(44),
        Tag::UmiQualityScores => Some(45),
        Tag::MateSequence => Some(46),
        Tag::ReadGroup => Some(47),
        Tag::UmiSequence => Some(48),
        Tag::OtherAlignments => Some(49),
        Tag::TemplateMappingQuality => Some(50),
        Tag::SegmentCount => Some(51),
        Tag::TranscriptStrand => Some(52),
        Tag::NextHitQualityScores => Some(53),
        Tag::SegmentLikelihood => Some(54),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let data = Data::try_from(vec![
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
            Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
        ])?;

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("".parse(), Ok(Data::default()));

        assert_eq!(
            "RG:Z:rg0\tNH:i:1".parse(),
            Ok(Data::try_from(vec![
                Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
                Field::new(Tag::AlignmentHitCount, Value::UInt8(1)),
            ])?)
        );

        assert_eq!(
            "NH:i:1\tNH:i:1".parse::<Data>(),
            Err(ParseError::DuplicateTag(Tag::AlignmentHitCount))
        );

        Ok(())
    }
}
