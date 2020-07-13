mod key;
mod ty;

pub use self::ty::Type;

use std::{convert::TryFrom, error, fmt};

use crate::record::genotype;

use super::{number, record, Number, Record};

use self::key::Key;

/// A VCF header genotype format record (`FORMAT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format {
    id: genotype::field::Key,
    number: Number,
    ty: Type,
    description: String,
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
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Key::Format.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, ",{}={}", Key::Number, self.number)?;
        write!(f, ",{}={}", Key::Type, self.ty)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a genotype format header
/// record.
#[derive(Debug)]
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
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Id))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Id) => v.parse().map_err(TryFromRecordError::InvalidId),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let number = it
        .next()
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Number))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Number) => v.parse().map_err(TryFromRecordError::InvalidNumber),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let ty = it
        .next()
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Type))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Type) => v.parse().map_err(TryFromRecordError::InvalidType),
            _ => Err(TryFromRecordError::MissingField(Key::Type)),
        })?;

    let description = it
        .next()
        .ok_or_else(|| TryFromRecordError::MissingField(Key::Description))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Description) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Description)),
        })?;

    Ok(Format {
        id,
        number,
        ty,
        description,
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
                (String::from("Type"), String::from("Integer")),
                (String::from("Description"), String::from("Genotype")),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let format = Format::try_from(record)?;

        let expected = r#"##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">"#;

        assert_eq!(format.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_format() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let format = Format::try_from(record)?;

        assert_eq!(format.id(), &genotype::field::Key::Genotype);
        assert_eq!(format.number(), Number::Count(1));
        assert_eq!(format.ty(), Type::Integer);
        assert_eq!(format.description(), "Genotype");

        Ok(())
    }
}
