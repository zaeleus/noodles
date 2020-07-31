//! SAM header reference sequence and fields.

mod molecule_topology;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt, num};

pub use self::{molecule_topology::MoleculeTopology, tag::Tag};

use super::{record, Record};

/// A SAM header reference sequence.
///
/// The reference sequence describes a sequence a read possibly mapped to. Both the reference
/// sequence name and length are guaranteed to be set.
///
/// A list of reference sequences creates a reference sequence dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    name: String,
    len: i32,
    fields: HashMap<Tag, String>,
}

#[allow(clippy::len_without_is_empty)]
impl ReferenceSequence {
    /// Creates a reference sequence with a name and length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// assert_eq!(reference_sequence.name(), "sq0");
    /// assert_eq!(reference_sequence.len(), 13);
    /// ```
    pub fn new(name: String, len: i32) -> Self {
        Self {
            name,
            len,
            fields: HashMap::new(),
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// assert_eq!(reference_sequence.name(), "sq0");
    /// ```
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// assert_eq!(reference_sequence.name(), "sq0");
    ///
    /// *reference_sequence.name_mut() = String::from("sq1");
    /// assert_eq!(reference_sequence.name(), "sq1");
    /// ```
    pub fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }

    /// Returns the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// assert_eq!(reference_sequence.len(), 13);
    /// ```
    pub fn len(&self) -> i32 {
        self.len
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// assert_eq!(reference_sequence.len(), 13);
    ///
    /// *reference_sequence.len_mut() = 8;
    /// assert_eq!(reference_sequence.len(), 8);
    /// ```
    pub fn len_mut(&mut self) -> &mut i32 {
        &mut self.len
    }

    /// Returns the raw fields of the reference sequence.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the name and length fields, as they are parsed and available as
    /// [`name`] and [`len`], respectively.
    ///
    /// [`name`]: #method.name
    /// [`len`]: #method.len
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Tag, ReferenceSequence};
    ///
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// reference_sequence.insert(Tag::Md5Checksum, String::from("d7eba311421bbc9d3ada44709dd61534"));
    ///
    /// let fields = reference_sequence.fields();
    ///
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(
    ///     fields.get(&Tag::Md5Checksum),
    ///     Some(&String::from("d7eba311421bbc9d3ada44709dd61534"))
    /// );
    ///
    /// assert_eq!(fields.get(&Tag::Name), None);
    /// assert_eq!(reference_sequence.name(), "sq0");
    ///
    /// assert_eq!(fields.get(&Tag::Length), None);
    /// assert_eq!(reference_sequence.len(), 13);
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }

    /// Returns a reference to the raw field value mapped to the given key.
    ///
    /// This can only be used for fields with unparsed values. For a reference sequence, [`name`]
    /// and [`len`] must be used instead of `get(Tag::Name)` and `get(Tag::Length)`, respectively.
    ///
    /// [`name`]: #method.name
    /// [`len`]: #method.len
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Tag, ReferenceSequence};
    ///
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// reference_sequence.insert(Tag::Md5Checksum, String::from("d7eba311421bbc9d3ada44709dd61534"));
    ///
    /// assert_eq!(
    ///     reference_sequence.get(&Tag::Md5Checksum),
    ///     Some(&String::from("d7eba311421bbc9d3ada44709dd61534"))
    /// );
    /// assert_eq!(reference_sequence.get(&Tag::AssemblyId), None);
    /// ```
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }

    /// Inserts a tag-raw value pair into the reference sequence.
    ///
    /// This follows similar semantics to [`std::collections::HashMap::insert`].
    ///
    /// [`std::collections::HashMap::insert`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.insert
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Tag, ReferenceSequence};
    /// let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);
    /// reference_sequence.insert(Tag::Md5Checksum, String::from("d7eba311421bbc9d3ada44709dd61534"));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: String) -> Option<String> {
        self.fields.insert(tag, value)
    }
}

impl fmt::Display for ReferenceSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::ReferenceSequence)?;
        write!(f, "\t{}:{}", Tag::Name, self.name)?;
        write!(f, "\t{}:{}", Tag::Length, self.len)?;

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header reference sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
    /// The length tag (`LN`) has a invalid value.
    InvalidLength(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
            Self::InvalidLength(e) => write!(f, "invalid reference sequence length: {}", e),
        }
    }
}

impl TryFrom<Record> for ReferenceSequence {
    type Error = ParseError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::ReferenceSequence, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(ParseError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Vec<(String, String)>) -> Result<ReferenceSequence, ParseError> {
    let mut name = None;
    let mut len = None;
    let mut fields = HashMap::new();

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

        match tag {
            Tag::Name => {
                name = Some(value);
            }
            Tag::Length => {
                len = value.parse().map(Some).map_err(ParseError::InvalidLength)?;
            }
            _ => {
                fields.insert(tag, value);
            }
        }
    }

    Ok(ReferenceSequence {
        name: name.ok_or_else(|| ParseError::MissingRequiredTag(Tag::Name))?,
        len: len.ok_or_else(|| ParseError::MissingRequiredTag(Tag::Length))?,
        fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let mut reference_sequence = ReferenceSequence::new(String::from("sq0"), 13);

        reference_sequence.fields.insert(
            Tag::Md5Checksum,
            String::from("d7eba311421bbc9d3ada44709dd61534"),
        );

        let actual = format!("{}", reference_sequence);
        let expected = "@SQ\tSN:sq0\tLN:13\tM5:d7eba311421bbc9d3ada44709dd61534";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(ParseError::InvalidRecord)
        );
    }

    #[test]
    fn test_from_str_with_missing_name() {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::Map(vec![
                (String::from("LN"), String::from("1")),
                (
                    String::from("M5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534"),
                ),
            ]),
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(ParseError::MissingRequiredTag(Tag::Name))
        );
    }

    #[test]
    fn test_from_str_with_missing_length() {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::Map(vec![
                (String::from("SN"), String::from("sq0")),
                (
                    String::from("M5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534"),
                ),
            ]),
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(ParseError::MissingRequiredTag(Tag::Length))
        );
    }

    #[test]
    fn test_from_str_with_missing_name_and_length() {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::Map(vec![(
                String::from("M5"),
                String::from("d7eba311421bbc9d3ada44709dd61534"),
            )]),
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(ParseError::MissingRequiredTag(Tag::Name))
        );
    }

    #[test]
    fn test_from_str_with_invalid_length() {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::Map(vec![(String::from("LN"), String::from("thirteen"))]),
        );

        assert!(matches!(
            ReferenceSequence::try_from(record),
            Err(ParseError::InvalidLength(_))
        ));
    }
}
