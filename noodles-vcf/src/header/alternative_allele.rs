//! VCF header symbolic alternate allele record and key.

mod key;

use std::{convert::TryFrom, error, fmt};

use crate::record::alternate_bases::allele::{symbol, Symbol};

use super::{record, Record};

use self::key::Key;

/// A VCF header symbolic alternate allele record (`ALT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeAllele {
    id: Symbol,
    description: String,
}

impl AlternativeAllele {
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
    ///     String::from("Deletion"),
    /// );
    /// ```
    pub fn new(id: Symbol, description: String) -> Self {
        Self { id, description }
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
    ///     String::from("Deletion"),
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
    ///     String::from("Deletion"),
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
        f.write_str("##")?;
        f.write_str(record::Key::AlternativeAllele.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

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
    MissingField(Key),
    /// The ID is invalid.
    InvalidId(symbol::ParseError),
}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing {} field", key),
            Self::InvalidId(e) => write!(f, "invalid ID: {}", e),
        }
    }
}

impl error::Error for TryFromRecordError {}

impl TryFrom<Record> for AlternativeAllele {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::AlternativeAllele, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<AlternativeAllele, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Id))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Id) => v.parse().map_err(TryFromRecordError::InvalidId),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let description = it
        .next()
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Description))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Description) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Description)),
        })?;

    Ok(AlternativeAllele { id, description })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::AlternativeAllele,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("DEL")),
                (String::from("Description"), String::from("Deletion")),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let alternative_allele = AlternativeAllele::try_from(record)?;

        let expected = r#"##ALT=<ID=DEL,Description="Deletion">"#;
        assert_eq!(alternative_allele.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_filter() {
        let record = build_record();

        assert_eq!(
            AlternativeAllele::try_from(record),
            Ok(AlternativeAllele {
                id: Symbol::StructuralVariant(symbol::StructuralVariant::from(
                    symbol::structural_variant::Type::Deletion
                )),
                description: String::from("Deletion"),
            })
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_record_key() {
        let record = Record::new(
            record::Key::FileFormat,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("DEL")),
                (String::from("Description"), String::from("Deletion")),
            ]),
        );

        assert_eq!(
            AlternativeAllele::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_record_value() {
        let record = Record::new(
            record::Key::AlternativeAllele,
            record::Value::String(String::from("VCFv4.3")),
        );

        assert_eq!(
            AlternativeAllele::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_a_missing_field() {
        let record = Record::new(
            record::Key::AlternativeAllele,
            record::Value::Struct(vec![(
                String::from("Description"),
                String::from("Deletion"),
            )]),
        );

        assert_eq!(
            AlternativeAllele::try_from(record),
            Err(TryFromRecordError::MissingField(Key::Id))
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_id() {
        let record = Record::new(
            record::Key::AlternativeAllele,
            record::Value::Struct(vec![
                (String::from("ID"), String::new()),
                (String::from("Description"), String::from("Deletion")),
            ]),
        );

        assert!(matches!(
            AlternativeAllele::try_from(record),
            Err(TryFromRecordError::InvalidId(_))
        ));
    }
}
