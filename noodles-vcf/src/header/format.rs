//! VCF header genotype format record and key.

pub mod key;
pub mod ty;

pub use self::{key::Key, ty::Type};

use std::{convert::TryFrom, error, fmt};

use indexmap::IndexMap;

use crate::record::genotype;

use super::{number, record, Number, Record};

/// A VCF header genotype format record (`FORMAT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format {
    id: genotype::field::Key,
    number: Number,
    ty: Type,
    description: String,
    fields: IndexMap<String, String>,
}

impl Format {
    /// Creates a VCF header genotype format record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    /// ```
    pub fn new(id: genotype::field::Key, number: Number, ty: Type, description: String) -> Self {
        Self {
            id,
            number,
            ty,
            description,
            fields: IndexMap::new(),
        }
    }

    /// Returns the genotype field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    ///
    /// assert_eq!(format.id(), &Key::Genotype);
    /// ```
    pub fn id(&self) -> &genotype::field::Key {
        &self.id
    }

    /// Returns the cardinality of the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    ///
    /// assert_eq!(format.number(), Number::Count(1));
    /// ```
    pub fn number(&self) -> Number {
        self.number
    }

    /// Returns the type of the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    ///
    /// assert_eq!(format.ty(), Type::String);
    /// ```
    pub fn ty(&self) -> Type {
        self.ty
    }

    /// Returns the description of the genotype field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    ///
    /// assert_eq!(format.description(), "Genotype");
    /// ```
    pub fn description(&self) -> &str {
        &self.description
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID`, `Number`, `Type`, and `Description`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    ///
    /// assert!(format.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Format.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, ",{}={}", Key::Number, self.number)?;
        write!(f, ",{}={}", Key::Type, self.ty)?;

        write!(f, ",{}=", Key::Description)?;
        super::fmt::write_escaped_string(f, self.description())?;

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a genotype format header
/// record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required field is missing.
    MissingField(Key),
    /// The ID is invalid.
    InvalidId(genotype::field::key::ParseError),
    /// The number is invalid.
    InvalidNumber(number::ParseError),
    /// The type is invalid.
    InvalidType(ty::ParseError),
    /// The number for the given ID does not match the number in the reserved definition.
    NumberMismatch(Number, Number),
    /// The type for the given ID does not match the type in the reserved definition.
    TypeMismatch(Type, Type),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing {} field", key),
            Self::InvalidId(e) => write!(f, "invalid ID: {}", e),
            Self::InvalidNumber(e) => write!(f, "invalid number: {}", e),
            Self::InvalidType(e) => write!(f, "invalid type: {}", e),
            Self::NumberMismatch(actual, expected) => {
                write!(f, "number mismatch: expected {}, got {}", expected, actual)
            }
            Self::TypeMismatch(actual, expected) => {
                write!(f, "type mismatch: expected {}, got {}", expected, actual)
            }
        }
    }
}

impl TryFrom<Record> for Format {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::Format, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Format, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Id))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Id) => v.parse().map_err(TryFromRecordError::InvalidId),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let number = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Number))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Number) => v.parse().map_err(TryFromRecordError::InvalidNumber),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let ty = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Type))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Type) => v.parse().map_err(TryFromRecordError::InvalidType),
            _ => Err(TryFromRecordError::MissingField(Key::Type)),
        })?;

    if !matches!(id, genotype::field::Key::Other(..)) {
        if id.number() != number {
            return Err(TryFromRecordError::NumberMismatch(number, id.number()));
        }

        if id.ty() != ty {
            return Err(TryFromRecordError::TypeMismatch(ty, id.ty()));
        }
    }

    let description = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Description))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Description) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Description)),
        })?;

    Ok(Format {
        id,
        number,
        ty,
        description,
        fields: it.collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let format = Format::try_from(record)?;

        let expected = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
        assert_eq!(format.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_format() {
        let record = build_record();

        assert_eq!(
            Format::try_from(record),
            Ok(Format::new(
                genotype::field::Key::Genotype,
                Number::Count(1),
                Type::String,
                String::from("Genotype"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_extra_fields() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
                (String::from("Comment"), String::from("noodles")),
            ]),
        );

        assert_eq!(
            Format::try_from(record),
            Ok(Format {
                id: genotype::field::Key::Genotype,
                number: Number::Count(1),
                ty: Type::String,
                description: String::from("Genotype"),
                fields: vec![(String::from("Comment"), String::from("noodles"))]
                    .into_iter()
                    .collect(),
            })
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_record_key() {
        let record = Record::new(
            record::Key::FileFormat,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert_eq!(
            Format::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_record_value() {
        let record = Record::new(
            record::Key::Format,
            record::Value::String(String::from("VCFv4.3")),
        );

        assert_eq!(
            Format::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_a_missing_field() {
        let record = Record::new(record::Key::Format, record::Value::Struct(Vec::new()));

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::MissingField(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_id() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::InvalidId(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_number() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("NA")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::InvalidNumber(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_type() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("str")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::InvalidType(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_a_number_mismatch() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from(".")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::NumberMismatch(
                Number::Unknown,
                Number::Count(1)
            ))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_a_type_mismatch() {
        let record = Record::new(
            record::Key::Format,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("GT")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        );

        assert!(matches!(
            Format::try_from(record),
            Err(TryFromRecordError::TypeMismatch(
                Type::Integer,
                Type::String
            ))
        ));
    }
}
