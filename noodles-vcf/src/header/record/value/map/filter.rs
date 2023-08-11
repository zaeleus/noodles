//! Inner VCF header filter map value.

pub(crate) mod tag;

pub use self::tag::Tag;

use std::{error, fmt, num};

use indexmap::IndexMap;

use self::tag::StandardTag;
use super::{builder, Described, Fields, Indexed, Inner, Map, OtherFields};

/// An inner VCF header filter map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filter {
    pub(crate) description: String,
    pub(crate) idx: Option<usize>,
}

impl Inner for Filter {
    type StandardTag = StandardTag;
    type Builder = builder::DescribedIndexed;
}

impl Described for Filter {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Indexed for Filter {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Filter> {
    /// Creates a default VCF header filter map value for "PASS".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Filter, Map};
    /// let actual = Map::<Filter>::pass();
    /// let expected = Map::<Filter>::new("All filters passed");
    /// assert_eq!(actual, expected);
    /// ```
    pub fn pass() -> Self {
        Self {
            inner: Filter {
                description: String::from("All filters passed"),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }

    /// Creates a VCF header filter map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Filter, Map};
    /// let map = Map::<Filter>::new("Quality below 10");
    /// ```
    pub fn new<D>(description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: Filter {
                description: description.into(),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }
}

impl fmt::Display for Map<Filter> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;

        if let Some(idx) = self.idx() {
            super::fmt_display_idx_field(f, idx)?;
        }

        Ok(())
    }
}

/// An error returned when a raw FILTER record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The IDX is invalid.
    InvalidIdx(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
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
            Self::InvalidIdx(_) => write!(f, "invalid IDX"),
        }
    }
}

impl TryFrom<Fields> for Map<Filter> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut description = None;
        let mut idx = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                tag::DESCRIPTION => try_replace(&mut description, tag::DESCRIPTION, value)?,
                tag::IDX => parse_idx(&value).and_then(|v| try_replace(&mut idx, tag::IDX, v))?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let description = description.ok_or(ParseError::MissingField(tag::DESCRIPTION))?;

        Ok(Self {
            inner: Filter { description, idx },
            other_fields,
        })
    }
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

impl builder::Inner<Filter> for builder::DescribedIndexed {
    fn build(self) -> Result<Filter, builder::BuildError> {
        let description = self
            .description
            .ok_or(builder::BuildError::MissingField("Description"))?;

        Ok(Filter {
            description,
            idx: self.idx,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let map = Map::<Filter>::pass();
        let expected = r#",Description="All filters passed""#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_filter() -> Result<(), ParseError> {
        let actual = Map::<Filter>::try_from(vec![(
            String::from("Description"),
            String::from("All filters passed"),
        )])?;

        let expected = Map::<Filter>::pass();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_filter_with_missing_fields() {
        assert_eq!(
            Map::<Filter>::try_from(vec![(String::from("Other"), String::from("noodles")),]),
            Err(ParseError::MissingField(tag::DESCRIPTION)),
        );
    }
}
