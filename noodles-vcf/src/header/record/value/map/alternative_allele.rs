mod builder;

use std::fmt;

use indexmap::IndexMap;

use super::{tag, Described, Fields, Inner, Map, TryFromFieldsError};
use crate::record::alternate_bases::allele::Symbol;

type StandardTag = tag::Described;
type Tag = tag::Tag<StandardTag>;

/// An inner VCF header alternative allele map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeAllele {
    description: String,
}

impl Inner for AlternativeAllele {
    type Id = Symbol;
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
    /// let map = Map::<AlternativeAllele>::new(
    ///     Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
    ///     "Deletion",
    /// );
    /// ```
    pub fn new<D>(id: Symbol, description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            id,
            inner: AlternativeAllele {
                description: description.into(),
            },
            other_fields: IndexMap::new(),
        }
    }
}

impl fmt::Display for Map<AlternativeAllele> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_prefix(f, self.id())?;
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;
        super::fmt_display_suffix(f)?;
        Ok(())
    }
}

impl TryFrom<Fields> for Map<AlternativeAllele> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields(fields.len());

        let mut id = None;
        let mut description = None;

        for (key, value) in fields {
            match Tag::from(key) {
                Tag::Standard(StandardTag::Id) => super::parse_id(&value, &mut id)?,
                Tag::Standard(StandardTag::Description) => {
                    super::parse_description(value, &mut description)?
                }
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let id = id.ok_or(TryFromFieldsError::MissingField("ID"))?;
        let description = description.ok_or(TryFromFieldsError::MissingField("Description"))?;

        Ok(Self {
            id,
            inner: AlternativeAllele { description },
            other_fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn del() -> Symbol {
        use crate::record::alternate_bases::allele::symbol::{
            structural_variant::Type, StructuralVariant,
        };

        Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion))
    }

    #[test]
    fn test_fmt() {
        let map = Map::<AlternativeAllele>::new(del(), "Deletion");
        let expected = r#"<ID=DEL,Description="Deletion">"#;
        assert_eq!(map.to_string(), expected);
    }

    #[test]
    fn test_try_from_fields_for_map_alternative_allele() -> Result<(), TryFromFieldsError> {
        let actual = Map::<AlternativeAllele>::try_from(vec![
            (String::from("ID"), String::from("DEL")),
            (String::from("Description"), String::from("Deletion")),
        ])?;

        let expected = Map::<AlternativeAllele>::new(del(), "Deletion");

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_alternative_allele_with_missing_fields() {
        assert_eq!(
            Map::<AlternativeAllele>::try_from(vec![(
                String::from("Description"),
                String::from("Deletion")
            ),]),
            Err(TryFromFieldsError::MissingField("ID")),
        );

        assert_eq!(
            Map::<AlternativeAllele>::try_from(vec![(String::from("ID"), String::from("DEL")),]),
            Err(TryFromFieldsError::MissingField("Description")),
        );
    }
}
