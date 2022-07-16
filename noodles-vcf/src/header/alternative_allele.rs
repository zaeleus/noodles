//! VCF header symbolic alternate allele record and key.

use std::{error, fmt};

use indexmap::IndexMap;

use super::record;
use crate::record::alternate_bases::allele::{symbol, Symbol};

const ID: &str = "ID";
const DESCRIPTION: &str = "Description";

/// A VCF header symbolic alternate allele record (`ALT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeAllele {
    id: Symbol,
    description: String,
}

impl AlternativeAllele {
    pub(crate) fn try_from_fields(
        id: String,
        fields: IndexMap<String, String>,
    ) -> Result<Self, TryFromRecordError> {
        parse_struct(id, fields)
    }

    /// Creates a VCF header symbolic alternate allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::AlternativeAllele,
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let alternative_allele = AlternativeAllele::new(
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
            description: description.into(),
        }
    }

    /// Returns the alternate allele symbol.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::AlternativeAllele,
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let alternative_allele = AlternativeAllele::new(
    ///     Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
    ///     "Deletion",
    /// );
    ///
    /// assert_eq!(
    ///     alternative_allele.id(),
    ///     &Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
    /// );
    /// ```
    pub fn id(&self) -> &Symbol {
        &self.id
    }

    /// Returns the description of the alternate allele symbol.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::AlternativeAllele,
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let alternative_allele = AlternativeAllele::new(
    ///     Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
    ///     "Deletion",
    /// );
    ///
    /// assert_eq!(alternative_allele.description(), "Deletion");
    /// ```
    pub fn description(&self) -> &str {
        &self.description
    }
}

impl fmt::Display for AlternativeAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::key::ALTERNATIVE_ALLELE.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id)?;

        write!(f, ",{}=", DESCRIPTION)?;
        super::fmt::write_escaped_string(f, self.description())?;

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to an alternative allele
/// header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required field is missing.
    MissingField(&'static str),
    /// The ID is invalid.
    InvalidId(symbol::ParseError),
}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing field: {}", key),
            Self::InvalidId(e) => write!(f, "invalid ID: {}", e),
        }
    }
}

impl error::Error for TryFromRecordError {}

fn parse_struct(
    raw_id: String,
    fields: IndexMap<String, String>,
) -> Result<AlternativeAllele, TryFromRecordError> {
    let id = raw_id.parse().map_err(TryFromRecordError::InvalidId)?;

    let mut it = fields.into_iter();

    let description = it
        .next()
        .ok_or(TryFromRecordError::MissingField(DESCRIPTION))
        .and_then(|(k, v)| match k.as_ref() {
            DESCRIPTION => Ok(v),
            _ => Err(TryFromRecordError::MissingField(DESCRIPTION)),
        })?;

    Ok(AlternativeAllele { id, description })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn del() -> Symbol {
        Symbol::StructuralVariant(symbol::StructuralVariant::from(
            symbol::structural_variant::Type::Deletion,
        ))
    }

    fn build_record() -> (String, IndexMap<String, String>) {
        (
            String::from("DEL"),
            [(String::from("Description"), String::from("Deletion"))]
                .into_iter()
                .collect(),
        )
    }

    #[test]
    fn test_fmt() {
        let alternative_allele = AlternativeAllele::new(del(), "Deletion");
        let expected = r#"##ALT=<ID=DEL,Description="Deletion">"#;
        assert_eq!(alternative_allele.to_string(), expected);
    }

    #[test]
    fn test_try_from_record_for_filter() {
        let (id, fields) = build_record();

        assert_eq!(
            AlternativeAllele::try_from_fields(id, fields),
            Ok(AlternativeAllele::new(del(), "Deletion"))
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_a_missing_field(
    ) -> Result<(), record::value::TryFromFieldsError> {
        let id = String::from("DEL");
        let fields = IndexMap::new();

        assert_eq!(
            AlternativeAllele::try_from_fields(id, fields),
            Err(TryFromRecordError::MissingField(DESCRIPTION))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_id(
    ) -> Result<(), record::value::TryFromFieldsError> {
        let id = String::new();
        let fields = [(String::from("Description"), String::from("Deletion"))]
            .into_iter()
            .collect();

        assert!(matches!(
            AlternativeAllele::try_from_fields(id, fields),
            Err(TryFromRecordError::InvalidId(_))
        ));

        Ok(())
    }
}
