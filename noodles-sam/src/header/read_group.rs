//! SAM header read group and fields.

mod platform;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{platform::Platform, tag::Tag};

use super::record;

/// A SAM header read group.
///
/// A read group typically defines the set of reads that came from the same run on a sequencing
/// instrument. The read group ID is guaranteed to be set.
#[derive(Debug)]
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
            ..Default::default()
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

impl Default for ReadGroup {
    fn default() -> Self {
        Self {
            id: String::new(),
            fields: HashMap::new(),
        }
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
#[derive(Debug)]
pub enum ParseError {
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<&[(String, String)]> for ReadGroup {
    type Error = ParseError;

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut read_group = ReadGroup::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

            if let Tag::Id = tag {
                read_group.id = value.into();
                has_id = true;
                continue;
            }

            read_group.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(ParseError::MissingRequiredTag(Tag::Id));
        }

        Ok(read_group)
    }
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
    fn test_from_str_with_no_version() {
        let fields = [(String::from("DS"), String::from("noodles"))];
        assert!(ReadGroup::try_from(&fields[..]).is_err());
    }
}
