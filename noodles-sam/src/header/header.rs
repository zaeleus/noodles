//! SAM header header and fields.
//!
//! The namespace of this module is intentionally awkward to disambiguate a SAM header
//! ([`crate::Header`]) and a header record ([`crate::header::header::Header`]).

mod builder;
mod group_order;
mod sort_order;
mod subsort_order;
mod tag;
mod version;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{
    builder::Builder, group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder,
    tag::Tag, version::Version,
};

use super::{record, Record};

/// A SAM header header.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    version: Version,
    sort_order: Option<SortOrder>,
    group_order: Option<GroupOrder>,
    subsort_order: Option<SubsortOrder>,
    fields: HashMap<Tag, String>,
}

impl Header {
    /// Creates a SAM header header builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let builder = Header::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a header with a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, Version};
    /// let header = Header::new(Version::new(1, 6));
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn new(version: Version) -> Self {
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
    /// use noodles_sam::header::header::{Header, Version};
    /// let header = Header::new(Version::new(1, 6));
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn version(&self) -> Version {
        self.version
    }

    /// Returns a mutable reference to the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, Version};
    ///
    /// let mut header = Header::new(Version::new(1, 6));
    /// assert_eq!(header.version(), Version::new(1, 6));
    ///
    /// *header.version_mut() = Version::new(1, 5);
    /// assert_eq!(header.version(), Version::new(1, 5));
    /// ```
    pub fn version_mut(&mut self) -> &mut Version {
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
    /// this will not include the version field, as it is parsed and available as
    /// [`Self::version`].
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, SortOrder, Tag, Version};
    ///
    /// let header = Header::builder()
    ///     .set_version(Version::new(1, 6))
    ///     .insert(Tag::Other(String::from("zn")), String::from("noodles"))
    ///     .build();
    ///
    /// let fields = header.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(
    ///     fields.get(&Tag::Other(String::from("zn"))),
    ///     Some(&String::from("noodles"))
    /// );
    ///
    /// assert_eq!(fields.get(&Tag::Version), None);
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }
}

impl Default for Header {
    fn default() -> Self {
        Builder::default().build()
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
    /// The version is invalid.
    InvalidVersion(version::ParseError),
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
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
            Self::InvalidVersion(e) => write!(f, "invalid version: {}", e),
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
    let mut builder = Header::builder();
    let mut version: Option<Version> = None;

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        builder = match tag {
            Tag::Version => {
                version = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidVersion)?;

                builder
            }
            Tag::SortOrder => {
                let sort_order = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidSortOrder)?;
                builder.set_sort_order(sort_order)
            }
            Tag::GroupOrder => {
                let group_order = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidGroupOrder)?;
                builder.set_group_order(group_order)
            }
            Tag::SubsortOrder => {
                let subsort_order = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidSubsortOrder)?;
                builder.set_subsort_order(subsort_order)
            }
            _ => builder.insert(tag, value),
        }
    }

    if let Some(v) = version {
        builder = builder.set_version(v);
    } else {
        return Err(TryFromRecordError::MissingRequiredTag(Tag::Version));
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Header::default();
        assert_eq!(header.version(), Version::default());
        assert!(header.sort_order().is_none());
        assert!(header.group_order().is_none());
        assert!(header.subsort_order().is_none());
        assert!(header.fields.is_empty());
    }

    #[test]
    fn test_fmt() {
        let header = Header::builder()
            .set_version(Version::new(1, 6))
            .set_sort_order(SortOrder::Unknown)
            .build();

        assert_eq!(header.to_string(), "@HD\tVN:1.6\tSO:unknown");
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
