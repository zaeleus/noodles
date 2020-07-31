//! SAM header read group and fields.

mod platform;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{platform::Platform, tag::Tag};

use super::{record, Record};

/// A SAM header read group.
///
/// A read group typically defines the set of reads that came from the same run on a sequencing
/// instrument. The read group ID is guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadGroup {
    id: String,
    fields: HashMap<Tag, String>,
}

impl ReadGroup {
    /// Creates a read group with an ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new(String::from("rg0"));
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn new(id: String) -> Self {
        Self {
            id,
            fields: HashMap::new(),
        }
    }

    /// Returns the read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new(String::from("rg0"));
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns a mutable reference to the read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let mut read_group = ReadGroup::new(String::from("rg0"));
    /// assert_eq!(read_group.id(), "rg0");
    ///
    /// *read_group.id_mut() = String::from("rg1");
    /// assert_eq!(read_group.id(), "rg1");
    /// ```
    //
    pub fn id_mut(&mut self) -> &mut String {
        &mut self.id
    }

    /// Returns the raw fields of the read group.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the ID field, as it is parsed and available as [`id`].
    ///
    /// [`id`]: #method.id
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{read_group::Tag, ReadGroup};
    ///
    /// let mut read_group = ReadGroup::new(String::from("rg0"));
    /// read_group.insert(Tag::Program, String::from("noodles-sam"));
    ///
    /// let fields = read_group.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(fields.get(&Tag::Program), Some(&String::from("noodles-sam")));
    ///
    /// assert_eq!(fields.get(&Tag::Id), None);
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }

    /// Returns a reference to the raw field value mapped to the given key.
    ///
    /// This can only be used for fields with unparsed values. For a read group, [`id`] must be
    /// used instead of `get(Tag::Id)`.
    ///
    /// [`id`]: #method.id
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{read_group::Tag, ReadGroup};
    ///
    /// let mut read_group = ReadGroup::new(String::from("rg0"));
    /// read_group.insert(Tag::Program, String::from("noodles-sam"));
    ///
    /// assert_eq!(read_group.get(&Tag::Program), Some(&String::from("noodles-sam")));
    /// assert_eq!(read_group.get(&Tag::Id), None);
    /// ```
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }

    /// Inserts a tag-raw value pair into the read group.
    ///
    /// This follows similar semantics to [`std::collections::HashMap::insert`].
    ///
    /// [`std::collections::HashMap::insert`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.insert
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{read_group::Tag, ReadGroup};
    /// let mut read_group = ReadGroup::new(String::from("rg0"));
    /// read_group.insert(Tag::Program, String::from("noodles-sam"));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: String) -> Option<String> {
        self.fields.insert(tag, value)
    }
}

impl fmt::Display for ReadGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::ReadGroup)?;
        write!(f, "\t{}:{}", Tag::Id, self.id)?;

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header read group fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<Record> for ReadGroup {
    type Error = ParseError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::ReadGroup, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(ParseError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Vec<(String, String)>) -> Result<ReadGroup, ParseError> {
    let mut id = None;
    let mut fields = HashMap::new();

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

        if let Tag::Id = tag {
            id = Some(value);
        } else {
            fields.insert(tag, value);
        }
    }

    Ok(ReadGroup {
        id: id.ok_or_else(|| ParseError::MissingRequiredTag(Tag::Id))?,
        fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let mut read_group = ReadGroup::new(String::from("rg0"));

        read_group
            .fields
            .insert(Tag::Program, String::from("noodles"));

        let actual = format!("{}", read_group);
        let expected = "@RG\tID:rg0\tPG:noodles";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(ReadGroup::try_from(record), Err(ParseError::InvalidRecord));
    }

    #[test]
    fn test_from_str_with_no_id() {
        let record = Record::new(
            record::Kind::ReadGroup,
            record::Value::Map(vec![(String::from("PG"), String::from("noodles"))]),
        );

        assert_eq!(
            ReadGroup::try_from(record),
            Err(ParseError::MissingRequiredTag(Tag::Id))
        );
    }
}
