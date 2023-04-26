//! Inner VCF header other map value.

pub(crate) mod tag;

pub use self::tag::Tag;

use std::{error, fmt};

use self::tag::StandardTag;
use super::{builder, Fields, Inner, Map, OtherFields};

/// An inner VCF header other map value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Other;

impl Inner for Other {
    type StandardTag = StandardTag;
    type Builder = builder::Identity;
}

impl Map<Other> {
    /// Creates a nonstandard VCF header map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Other, Map};
    /// let map = Map::<Other>::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }
}

impl fmt::Display for Map<Other> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_other_fields(f, self.other_fields())
    }
}

/// An error returned when a raw other record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

impl TryFrom<Fields> for Map<Other> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        Ok(Self {
            inner: Other,
            other_fields,
        })
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
    fn test_fmt() {
        let map = Map::<Other>::new();
        assert!(map.to_string().is_empty());
    }

    #[test]
    fn test_try_from_fields_for_map_other() -> Result<(), Box<dyn std::error::Error>> {
        let actual = Map::<Other>::try_from(vec![(String::from("noodles"), String::from("vcf"))])?;

        let expected = Map::<Other>::builder()
            .insert("noodles".parse()?, "vcf")
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
