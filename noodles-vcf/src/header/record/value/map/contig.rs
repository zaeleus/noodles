//! Inner VCF header contig map value.

mod builder;
pub mod name;
mod tag;

pub use self::{name::Name, tag::Tag};

use std::fmt;

use self::tag::StandardTag;
use super::{Fields, Indexed, Inner, Map, OtherFields, TryFromFieldsError};

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

impl TryFrom<Fields> for Map<Contig> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut length = None;
        let mut md5 = None;
        let mut url = None;
        let mut idx = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
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

fn parse_length(s: &str) -> Result<usize, TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("length"))
}

fn parse_idx(s: &str) -> Result<usize, TryFromFieldsError> {
    s.parse()
        .map_err(|_| TryFromFieldsError::InvalidValue("IDX"))
}

fn try_replace<T>(option: &mut Option<T>, _: Tag, value: T) -> Result<(), TryFromFieldsError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(TryFromFieldsError::DuplicateTag)
    }
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), TryFromFieldsError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(_) => Err(TryFromFieldsError::DuplicateTag),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), TryFromFieldsError> {
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
