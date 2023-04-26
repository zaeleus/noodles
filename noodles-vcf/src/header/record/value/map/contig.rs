//! Inner VCF header contig map value.

mod builder;
pub mod name;
pub(crate) mod tag;

pub use self::{name::Name, tag::Tag};

use std::{error, fmt, num};

use self::tag::StandardTag;
use super::{Fields, Indexed, Inner, Map, OtherFields};

/// An inner VCF header contig map value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Contig {
    length: Option<usize>,
    md5: Option<String>,
    url: Option<String>,
    idx: Option<usize>,
}

impl Inner for Contig {
    type StandardTag = StandardTag;
    type Builder = builder::Builder;
}

impl Indexed for Contig {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Contig> {
    /// Creates a VCF header contig map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.length().is_none());
    /// ```
    pub fn length(&self) -> Option<usize> {
        self.inner.length
    }

    /// Returns a mutable reference to the length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.length().is_none());
    ///
    /// *map.length_mut() = Some(8);
    /// assert_eq!(map.length(), Some(8));
    /// ```
    pub fn length_mut(&mut self) -> &mut Option<usize> {
        &mut self.inner.length
    }

    /// Returns the MD5 hexdigest.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.md5().is_none());
    /// ```
    pub fn md5(&self) -> Option<&str> {
        self.inner.md5.as_deref()
    }

    /// Returns a mutable reference to the MD5 hexdigest.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.md5().is_none());
    ///
    /// *map.md5_mut() = Some(String::from("d7eba311421bbc9d3ada44709dd61534"));
    /// assert_eq!(map.md5(), Some("d7eba311421bbc9d3ada44709dd61534"));
    /// ```
    pub fn md5_mut(&mut self) -> &mut Option<String> {
        &mut self.inner.md5
    }

    /// Returns the URL.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.url().is_none());
    /// ```
    pub fn url(&self) -> Option<&str> {
        self.inner.url.as_deref()
    }

    /// Returns a mutable reference to the URL.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.url().is_none());
    ///
    /// *map.url_mut() = Some(String::from("https://example.com/reference.fa"));
    /// assert_eq!(map.url(), Some("https://example.com/reference.fa"));
    /// ```
    pub fn url_mut(&mut self) -> &mut Option<String> {
        &mut self.inner.url
    }
}

impl fmt::Display for Map<Contig> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(length) = self.length() {
            write!(f, ",{tag}={length}", tag = tag::LENGTH)?;
        }

        if let Some(md5) = self.md5() {
            write!(f, ",{tag}={md5}", tag = tag::MD5)?;
        }

        if let Some(url) = self.url() {
            write!(f, ",{tag}={url}", tag = tag::URL)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        if let Some(idx) = self.idx() {
            super::fmt_display_idx_field(f, idx)?;
        }

        Ok(())
    }
}

/// An error returned when a raw contig record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The ID is invalid.
    InvalidId(name::ParseError),
    /// The length is invalid.
    InvalidLength(num::ParseIntError),
    /// The IDX is invalid.
    InvalidIdx(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
            Self::InvalidIdx(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::InvalidIdx(_) => write!(f, "invalid IDX"),
        }
    }
}

impl TryFrom<Fields> for Map<Contig> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut length = None;
        let mut md5 = None;
        let mut url = None;
        let mut idx = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                tag::LENGTH => {
                    parse_length(&value).and_then(|v| try_replace(&mut length, tag::LENGTH, v))?
                }
                tag::MD5 => try_replace(&mut md5, tag::MD5, value)?,
                tag::URL => try_replace(&mut url, tag::URL, value)?,
                tag::IDX => parse_idx(&value).and_then(|v| try_replace(&mut idx, tag::IDX, v))?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        Ok(Self {
            inner: Contig {
                length,
                md5,
                url,
                idx,
            },
            other_fields,
        })
    }
}

fn parse_length(s: &str) -> Result<usize, ParseError> {
    s.parse().map_err(ParseError::InvalidLength)
}

fn parse_idx(s: &str) -> Result<usize, ParseError> {
    s.parse().map_err(ParseError::InvalidIdx)
}

fn try_replace<T>(option: &mut Option<T>, tag: Tag, value: T) -> Result<(), ParseError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(ParseError::DuplicateTag(tag))
    }
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => {
            let (t, _) = entry.remove_entry();
            Err(ParseError::DuplicateTag(Tag::Other(t)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let map = Map::<Contig>::try_from(vec![
            (String::from("length"), String::from("8")),
            (
                String::from("md5"),
                String::from("d7eba311421bbc9d3ada44709dd61534"),
            ),
            (
                String::from("URL"),
                String::from("https://example.com/reference.fa"),
            ),
        ])?;

        let expected = r#",length=8,md5=d7eba311421bbc9d3ada44709dd61534,URL=https://example.com/reference.fa"#;
        assert_eq!(map.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_contig() -> Result<(), Box<dyn std::error::Error>> {
        let actual = Map::<Contig>::try_from(Vec::new())?;
        let expected = Map::<Contig>::new();
        assert_eq!(actual, expected);
        Ok(())
    }
}
