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

use super::record;

static VERSION: &str = "1.6";

/// A SAM header header.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    version: String,
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
            fields: HashMap::new(),
        }
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::Header)?;
        write!(f, "\t{}:{}", Tag::Version, self.version)?;

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
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

impl TryFrom<&[(String, String)]> for Header {
    type Error = ParseError;

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut header = Header::default();

        let mut has_version = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

            if let Tag::Version = tag {
                header.version = value.into();
                has_version = true;
                continue;
            }

            header.fields.insert(tag, value.into());
        }

        if !has_version {
            return Err(ParseError::MissingRequiredTag(Tag::Version));
        }

        Ok(header)
    }
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
    fn test_from_str_with_no_version() {
        let fields = [(String::from("SO"), String::from("coordinate"))];

        assert_eq!(
            Header::try_from(&fields[..]),
            Err(ParseError::MissingRequiredTag(Tag::Version))
        );
    }
}
