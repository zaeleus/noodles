//! VCF header genotype format record and key.

mod builder;
pub mod key;
mod tag;
pub mod ty;

pub use self::{key::Key, tag::Tag, ty::Type};

use std::{error, fmt, num};

use indexmap::IndexMap;

use self::builder::Builder;
use super::{number, record, FileFormat, Number};

/// A VCF header genotype format record (`FORMAT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format {
    id: Key,
    number: Number,
    ty: Type,
    description: String,
    idx: Option<usize>,
    fields: IndexMap<tag::Other, String>,
}

impl Format {
    pub(crate) fn try_from_fields(
        id: String,
        fields: IndexMap<String, String>,
        file_format: FileFormat,
    ) -> Result<Self, TryFromRecordError> {
        parse_struct(file_format, id, fields)
    }

    /// Creates a VCF header genotype format record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{format::{Key, Type}, Format, Number};
    ///
    /// let format = Format::new(
    ///     Key::Genotype,
    ///     Number::Count(1),
    ///     Type::String,
    ///     String::from("Genotype"),
    /// );
    /// ```
    pub fn new(id: Key, number: Number, ty: Type, description: String) -> Self {
        Self {
            id,
            number,
            ty,
            description,
            idx: None,
            fields: IndexMap::new(),
        }
    }

    /// Returns the genotype field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{format::Key, Format};
    /// let format = Format::from(Key::Genotype);
    /// assert_eq!(format.id(), &Key::Genotype);
    /// ```
    pub fn id(&self) -> &Key {
        &self.id
    }

    /// Returns the cardinality of the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{format::Key, Format, Number};
    /// let format = Format::from(Key::Genotype);
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
    /// use noodles_vcf::header::{format::{Key, Type}, Format};
    /// let format = Format::from(Key::Genotype);
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
    /// use noodles_vcf::header::{format::Key, Format};
    /// let format = Format::from(Key::Genotype);
    /// assert_eq!(format.description(), "Genotype");
    /// ```
    pub fn description(&self) -> &str {
        &self.description
    }

    /// Returns the index of the ID in the dictionary of strings.
    ///
    /// This is typically used in BCF.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{format::Key, Format};
    /// let format = Format::from(Key::Genotype);
    /// assert!(format.idx().is_none());
    /// ```
    pub fn idx(&self) -> Option<usize> {
        self.idx
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID`, `Number`, `Type`, `Description`, and `IDX`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{format::Key, Format};
    /// let format = Format::from(Key::Genotype);
    /// assert!(format.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<tag::Other, String> {
        &self.fields
    }
}

impl From<Key> for Format {
    fn from(key: Key) -> Self {
        let number = key::number(&key).unwrap_or(Number::Count(1));
        let ty = key::ty(&key).unwrap_or(Type::String);
        let description = key::description(&key).map(|s| s.into()).unwrap_or_default();
        Self::new(key, number, ty, description)
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::key::FORMAT.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", tag::ID, self.id)?;
        write!(f, ",{}={}", tag::NUMBER, self.number)?;
        write!(f, ",{}={}", tag::TYPE, self.ty)?;

        write!(f, ",{}=", tag::DESCRIPTION)?;
        super::fmt::write_escaped_string(f, self.description())?;

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        if let Some(idx) = self.idx() {
            write!(f, ",{}={}", tag::IDX, idx)?;
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
    /// The ID is invalid.
    InvalidId(key::ParseError),
    /// The number is invalid.
    InvalidNumber(number::ParseError),
    /// The type is invalid.
    InvalidType(ty::ParseError),
    /// The index (`IDX`) is invalid.
    InvalidIdx(num::ParseIntError),
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
            Self::InvalidId(e) => write!(f, "invalid ID: {}", e),
            Self::InvalidNumber(e) => write!(f, "invalid number: {}", e),
            Self::InvalidType(e) => write!(f, "invalid type: {}", e),
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", tag::IDX, e),
            Self::NumberMismatch(actual, expected) => {
                write!(f, "number mismatch: expected {}, got {}", expected, actual)
            }
            Self::TypeMismatch(actual, expected) => {
                write!(f, "type mismatch: expected {}, got {}", expected, actual)
            }
        }
    }
}

fn parse_struct(
    file_format: FileFormat,
    raw_id: String,
    fields: IndexMap<String, String>,
) -> Result<Format, TryFromRecordError> {
    let mut builder = Builder::default();

    let id = raw_id.parse().map_err(TryFromRecordError::InvalidId)?;
    builder = builder.set_id(id);

    for (key, value) in fields {
        builder = match Tag::from(key) {
            tag::ID => todo!(),
            tag::NUMBER => {
                let number = value.parse().map_err(TryFromRecordError::InvalidNumber)?;
                builder.set_number(number)
            }
            tag::TYPE => {
                let ty = value.parse().map_err(TryFromRecordError::InvalidType)?;
                builder.set_type(ty)
            }
            tag::DESCRIPTION => builder.set_description(value),
            tag::IDX => {
                let idx = value.parse().map_err(TryFromRecordError::InvalidIdx)?;
                builder.set_idx(idx)
            }
            Tag::Other(t) => builder.insert(t, value),
        };
    }

    let format = builder
        .build()
        .map_err(|_| TryFromRecordError::InvalidRecord)?;

    let id = format.id();

    if file_format >= FileFormat::new(4, 3) && !matches!(id, Key::Other(..)) {
        if let (Some(expected_number), Some(expected_type)) = (key::number(id), key::ty(id)) {
            if format.number() != expected_number {
                return Err(TryFromRecordError::NumberMismatch(
                    format.number(),
                    expected_number,
                ));
            }

            if format.ty() != expected_type {
                return Err(TryFromRecordError::TypeMismatch(format.ty(), expected_type));
            }
        }
    }

    Ok(format)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> (String, IndexMap<String, String>) {
        (
            String::from("GT"),
            [
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (String::from("Description"), String::from("Genotype")),
            ]
            .into_iter()
            .collect(),
        )
    }

    #[test]
    fn test_from_genotype_field_key_for_format() {
        let actual = Format::from(Key::Genotype);

        let expected = Format::new(
            Key::Genotype,
            Number::Count(1),
            Type::String,
            String::from("Genotype"),
        );

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let (id, fields) = build_record();
        let format = Format::try_from_fields(id, fields, Default::default())?;

        let expected = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
        assert_eq!(format.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_file_format() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from(".")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert_eq!(
            Format::try_from_fields(id, fields, FileFormat::new(4, 2)),
            Ok(Format::new(
                Key::Genotype,
                Number::Unknown,
                Type::String,
                String::from("Genotype"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_format() {
        let (id, fields) = build_record();

        assert_eq!(
            Format::try_from_fields(id, fields, Default::default()),
            Ok(Format::new(
                Key::Genotype,
                Number::Count(1),
                Type::String,
                String::from("Genotype"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_extra_fields() -> Result<(), &'static str> {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
            (String::from("Comment"), String::from("noodles")),
            (String::from("IDX"), String::from("1")),
        ]
        .into_iter()
        .collect();

        assert_eq!(
            Format::try_from_fields(id, fields, Default::default()),
            Ok(Format {
                id: Key::Genotype,
                number: Number::Count(1),
                ty: Type::String,
                description: String::from("Genotype"),
                idx: Some(1),
                fields: [(
                    Tag::other("Comment").ok_or("invalid tag")?,
                    String::from("noodles")
                )]
                .into_iter()
                .collect(),
            })
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_format_with_a_missing_field() {
        assert_eq!(
            Format::try_from_fields(String::from("GT"), IndexMap::new(), Default::default()),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_id() {
        let id = String::new();
        let fields = [
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::InvalidId(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_number() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from("one")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::InvalidNumber(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_type() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("str")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::InvalidType(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_an_invalid_idx() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
            (String::from("IDX"), String::from("zero")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_a_number_mismatch() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from(".")),
            (String::from("Type"), String::from("String")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::NumberMismatch(
                Number::Unknown,
                Number::Count(1)
            ))
        ));
    }

    #[test]
    fn test_try_from_record_for_format_with_a_type_mismatch() {
        let id = String::from("GT");
        let fields = [
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (String::from("Description"), String::from("Genotype")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Format::try_from_fields(id, fields, Default::default()),
            Err(TryFromRecordError::TypeMismatch(
                Type::Integer,
                Type::String
            ))
        ));
    }
}
