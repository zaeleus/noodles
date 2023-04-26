//! Inner VCF header alternative allele map value.

mod builder;
pub(crate) mod tag;

pub use self::tag::Tag;

use std::{error, fmt};

use self::tag::StandardTag;
use super::{Described, Fields, Inner, Map, OtherFields};

/// An inner VCF header alternative allele map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeAllele {
    description: String,
}

impl Inner for AlternativeAllele {
    type StandardTag = StandardTag;
    type Builder = builder::Builder;
}

impl Described for AlternativeAllele {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Map<AlternativeAllele> {
    /// Creates a VCF header alternative allele map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::record::value::{map::AlternativeAllele, Map},
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let id = Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion));
    /// let map = Map::<AlternativeAllele>::new("Deletion");
    /// ```
    pub fn new<D>(description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: AlternativeAllele {
                description: description.into(),
            },
            other_fields: OtherFields::new(),
        }
    }
}

impl fmt::Display for Map<AlternativeAllele> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;
        Ok(())
    }
}

/// An error returned when a raw ALT record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The ID is invalid.
    InvalidId(crate::record::alternate_bases::allele::symbol::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
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
        }
    }
}

impl TryFrom<Fields> for Map<AlternativeAllele> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut description = None;
        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(ParseError::DuplicateTag(tag::ID)),
                tag::DESCRIPTION => try_replace(&mut description, tag::DESCRIPTION, value)?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let description = description.ok_or(ParseError::MissingField(tag::DESCRIPTION))?;

        Ok(Self {
            inner: AlternativeAllele { description },
            other_fields,
        })
    }
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
    fn test_fmt() {
        let map = Map::<AlternativeAllele>::new("Deletion");
        let expected = r#",Description="Deletion""#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_alternative_allele() -> Result<(), ParseError> {
        let actual = Map::<AlternativeAllele>::try_from(vec![(
            String::from("Description"),
            String::from("Deletion"),
        )])?;

        let expected = Map::<AlternativeAllele>::new("Deletion");

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_alternative_allele_with_missing_fields() {
        assert_eq!(
            Map::<AlternativeAllele>::try_from(vec![(
                String::from("Other"),
                String::from("noodles")
            ),]),
            Err(ParseError::MissingField(tag::DESCRIPTION)),
        );
    }
}
