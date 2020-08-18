//! SAM header header and fields.
//!
//! The namespace of this module is intentionally awkward to disambiguate a SAM header
//! ([`sam::Header`]) and a header record ([`sam::header::header::Header`]).
//!
//! [`sam::Header`]: ../struct.Header.html
//! [`sam::header::header::Header`]: struct.Header.html

mod group_order;
mod sort_order;
mod subsort_order;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, tag::Tag,
};

use super::{record, Record};

static VERSION: &str = "1.6";

/// A SAM header header.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    version: String,
    sort_order: Option<SortOrder>,
    group_order: Option<GroupOrder>,
    subsort_order: Option<SubsortOrder>,
    fields: HashMap<Tag, String>,
}

impl Header {
    /// Creates a header with a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::new(String::from("1.6"));
    /// assert_eq!(header.version(), "1.6");
    /// ```
    pub fn new(version: String) -> Self {
        Self {
            version,
            ..Default::default()
        }
    }

    /// Returns the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::new(String::from("1.6"));
    /// assert_eq!(header.version(), "1.6");
    /// ```
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Returns a mutable reference to the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    ///
    /// let mut header = Header::new(String::from("1.6"));
    /// assert_eq!(header.version(), "1.6");
    ///
    /// *header.version_mut() = String::from("1.5");
    /// assert_eq!(header.version(), "1.5");
    /// ```
    pub fn version_mut(&mut self) -> &mut String {
        &mut self.version
    }

    /// Returns the sort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::default();
    /// assert!(header.sort_order().is_none());
    /// ```
    pub fn sort_order(&self) -> Option<SortOrder> {
        self.sort_order
    }

    /// Returns the group order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::default();
    /// assert!(header.group_order().is_none());
    /// ```
    pub fn group_order(&self) -> Option<GroupOrder> {
        self.group_order
    }

    /// Returns the subsort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::default();
    /// assert!(header.subsort_order().is_none());
    /// ```
    pub fn subsort_order(&self) -> Option<&SubsortOrder> {
        self.subsort_order.as_ref()
    }

    /// Returns the raw fields of the header.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the version field, as it is parsed and available as [`version`].
    ///
    /// [`version`]: #method.version
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{self, Header};
    ///
    /// let mut header = Header::new(String::from("1.6"));
    /// header.insert(header::Tag::SortOrder, String::from("coordinate"));
    ///
    /// let fields = header.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(fields.get(&header::Tag::SortOrder), Some(&String::from("coordinate")));
    /// assert_eq!(fields.get(&header::Tag::Version), None);
    /// assert_eq!(header.version(), "1.6");
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }

    /// Returns a reference to the raw field value mapped to the given key.
    ///
    /// This can only be used for fields with unparsed values. For the header, [`version`] must be
    /// used instead of `get(header::Tag::Version)`.
    ///
    /// [`version`]: #method.version
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{self, Header};
    ///
    /// let mut header = Header::default();
    /// header.insert(header::Tag::SortOrder, String::from("coordinate"));
    ///
    /// assert_eq!(header.get(&header::Tag::SortOrder), Some(&String::from("coordinate")));
    /// assert_eq!(header.get(&header::Tag::GroupOrder), None);
    /// ```
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }

    /// Inserts a tag-raw value pair into the header.
    ///
    /// This follows similar semantics to [`std::collections::HashMap::insert`].
    ///
    /// [`std::collections::HashMap::insert`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.insert
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{self, Header};
    /// let mut header = Header::default();
    /// header.insert(header::Tag::SortOrder, String::from("coordinate"));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: String) -> Option<String> {
        self.fields.insert(tag, value)
    }
}

impl Default for Header {
    fn default() -> Self {
        Header {
            version: VERSION.into(),
            sort_order: None,
            group_order: None,
            subsort_order: None,
            fields: HashMap::new(),
        }
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::Header)?;
        write!(f, "\t{}:{}", Tag::Version, self.version)?;

        if let Some(sort_order) = self.sort_order {
            write!(f, "\t{}:{}", Tag::SortOrder, sort_order)?;
        }

        if let Some(group_order) = self.group_order {
            write!(f, "\t{}:{}", Tag::GroupOrder, group_order)?;
        }

        if let Some(subsort_order) = self.group_order {
            write!(f, "\t{}:{}", Tag::SubsortOrder, subsort_order)?;
        }

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
    /// The sort order is invalid.
    InvalidSortOrder(sort_order::ParseError),
    /// The group order is invalid.
    InvalidGroupOrder(group_order::ParseError),
    /// The subsort order is invalid.
    InvalidSubsortOrder(subsort_order::ParseError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
            Self::InvalidSortOrder(e) => write!(f, "invalid sort order: {}", e),
            Self::InvalidGroupOrder(e) => write!(f, "invalid group order: {}", e),
            Self::InvalidSubsortOrder(e) => write!(f, "invalid subsort order: {}", e),
        }
    }
}

impl TryFrom<Record> for Header {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::Header, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Vec<(String, String)>) -> Result<Header, TryFromRecordError> {
    let mut version = None;
    let mut sort_order = None;
    let mut group_order = None;
    let mut subsort_order = None;
    let mut fields = HashMap::new();

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        match tag {
            Tag::Version => version = Some(value),
            Tag::SortOrder => {
                sort_order = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidSortOrder)?
            }
            Tag::GroupOrder => {
                group_order = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidGroupOrder)?
            }
            Tag::SubsortOrder => {
                subsort_order = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidSubsortOrder)?
            }
            _ => {
                fields.insert(tag, value);
            }
        }
    }

    Ok(Header {
        version: version.ok_or_else(|| TryFromRecordError::MissingRequiredTag(Tag::Version))?,
        sort_order,
        group_order,
        subsort_order,
        fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Header::default();
        assert_eq!(header.version(), "1.6");
        assert!(header.fields.is_empty());
    }

    #[test]
    fn test_fmt() {
        let mut header = Header::new(String::from("1.6"));

        header
            .fields
            .insert(Tag::SortOrder, String::from("unknown"));

        let actual = format!("{}", header);
        let expected = "@HD\tVN:1.6\tSO:unknown";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_try_from_record_for_header_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(
            Header::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_header_with_no_version() {
        let record = Record::new(
            record::Kind::Header,
            record::Value::Map(vec![(String::from("SO"), String::from("coordinate"))]),
        );

        assert_eq!(
            Header::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Version))
        );
    }
}
