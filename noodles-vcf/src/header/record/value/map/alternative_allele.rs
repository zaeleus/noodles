//! Inner VCF header alternative allele map value.

mod builder;
mod tag;

pub use self::tag::Tag;

use std::fmt;

use self::tag::StandardTag;
use super::{Described, Fields, Inner, Map, OtherFields, TryFromFieldsError};

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

impl TryFrom<Fields> for Map<AlternativeAllele> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut description = None;
        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
                tag::DESCRIPTION => try_replace(&mut description, tag::DESCRIPTION, value)?,
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let description = description.ok_or(TryFromFieldsError::MissingField("Description"))?;

        Ok(Self {
            inner: AlternativeAllele { description },
            other_fields,
        })
    }
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
    fn test_fmt() {
        let map = Map::<AlternativeAllele>::new("Deletion");
        let expected = r#",Description="Deletion""#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_alternative_allele() -> Result<(), TryFromFieldsError> {
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
            Err(TryFromFieldsError::MissingField("Description")),
        );
    }
}
