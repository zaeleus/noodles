//! VCF header information record and components.

mod builder;
pub mod key;
mod tag;
pub mod ty;

pub use self::{key::Key, tag::Tag, ty::Type};

use std::{error, fmt, num};

use indexmap::IndexMap;

use self::builder::Builder;
use super::{number, record, FileFormat, Number, Record};

const ID: &str = "ID";
const NUMBER: &str = "Number";
const TYPE: &str = "Type";
const DESCRIPTION: &str = "Description";
const IDX: &str = "IDX";

/// A VCF header information record (`INFO`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Info {
    id: Key,
    number: Number,
    ty: Type,
    description: String,
    idx: Option<usize>,
    fields: IndexMap<tag::Other, String>,
}

impl Info {
    /// Converts a generic VCF header record to a VCF header info record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{
    ///     info::{Key as InfoKey, Type},
    ///     record::{Key as RecordKey, Value},
    ///     FileFormat, Info, Number, Record,
    /// };
    ///
    /// let record = Record::new(
    ///     RecordKey::Info,
    ///     Value::Struct(vec![
    ///         (String::from("ID"), String::from("NS")),
    ///         (String::from("Number"), String::from("1")),
    ///         (String::from("Type"), String::from("Integer")),
    ///         (
    ///             String::from("Description"),
    ///             String::from("Number of samples with data"),
    ///         ),
    ///     ]),
    /// );
    ///
    /// assert_eq!(
    ///     Info::try_from_record_file_format(record, FileFormat::new(4, 3)),
    ///     Ok(Info::new(
    ///         InfoKey::SamplesWithDataCount,
    ///         Number::Count(1),
    ///         Type::Integer,
    ///         String::from("Number of samples with data"),
    ///     ))
    /// );
    /// ```
    pub fn try_from_record_file_format(
        record: Record,
        file_format: FileFormat,
    ) -> Result<Self, TryFromRecordError> {
        match record.into() {
            (record::Key::Info, record::Value::Struct(fields)) => parse_struct(file_format, fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }

    /// Creates a VCF header information record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{info::{Key, Type}, Info, Number};
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
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

    /// Returns the information field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{info::Key, Info};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert_eq!(info.id(), &Key::SamplesWithDataCount);
    /// ```
    pub fn id(&self) -> &Key {
        &self.id
    }

    /// Returns the cardinality of the information field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{info::Key, Info, Number};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert_eq!(info.number(), Number::Count(1));
    /// ```
    pub fn number(&self) -> Number {
        self.number
    }

    /// Returns the type of the information field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{info::{Key, Type}, Info};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert_eq!(info.ty(), Type::Integer);
    /// ```
    pub fn ty(&self) -> Type {
        self.ty
    }

    /// Returns the description of the information field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{info::Key, Info};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert_eq!(info.description(), "Number of samples with data");
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
    /// use noodles_vcf::header::{info::Key, Info};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert!(info.idx().is_none());
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
    /// use noodles_vcf::header::{info::Key, Info};
    /// let info = Info::from(Key::SamplesWithDataCount);
    /// assert!(info.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<tag::Other, String> {
        &self.fields
    }
}

impl From<Key> for Info {
    fn from(key: Key) -> Self {
        let number = key::number(&key).unwrap_or(Number::Count(1));
        let ty = key::ty(&key).unwrap_or(Type::String);
        let description = key::description(&key).map(|s| s.into()).unwrap_or_default();
        Self::new(key, number, ty, description)
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Info.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id)?;
        write!(f, ",{}={}", NUMBER, self.number)?;
        write!(f, ",{}={}", TYPE, self.ty)?;

        write!(f, ",{}=", DESCRIPTION)?;
        super::fmt::write_escaped_string(f, self.description())?;

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        if let Some(idx) = self.idx() {
            write!(f, ",{}={}", IDX, idx)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a info header record.
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
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", IDX, e),
            Self::NumberMismatch(actual, expected) => {
                write!(f, "number mismatch: expected {}, got {}", expected, actual)
            }
            Self::TypeMismatch(actual, expected) => {
                write!(f, "type mismatch: expected {}, got {}", expected, actual)
            }
        }
    }
}

impl TryFrom<Record> for Info {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        Self::try_from_record_file_format(record, FileFormat::default())
    }
}

fn parse_struct(
    file_format: FileFormat,
    fields: Vec<(String, String)>,
) -> Result<Info, TryFromRecordError> {
    let mut builder = Builder::default();

    for (key, value) in fields {
        let tag = Tag::from(key);

        builder = match tag {
            tag::ID => {
                let id = value.parse().map_err(TryFromRecordError::InvalidId)?;
                builder.set_id(id)
            }
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
        }
    }

    let info = builder
        .build()
        .map_err(|_| TryFromRecordError::InvalidRecord)?;

    let id = info.id();

    if file_format >= FileFormat::new(4, 3) && !matches!(id, Key::Other(..)) {
        if let (Some(expected_number), Some(expected_type)) = (key::number(id), key::ty(id)) {
            if info.number() != expected_number {
                return Err(TryFromRecordError::NumberMismatch(
                    info.number(),
                    expected_number,
                ));
            }

            if info.ty() != expected_type {
                return Err(TryFromRecordError::TypeMismatch(info.ty(), expected_type));
            }
        }
    }

    Ok(info)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        )
    }

    #[test]
    fn test_from_info_field_key_for_info() {
        let actual = Info::from(Key::SamplesWithDataCount);

        let expected = Info::new(
            Key::SamplesWithDataCount,
            Number::Count(1),
            Type::Integer,
            String::from("Number of samples with data"),
        );

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let info = Info::try_from(record)?;

        let expected =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;
        assert_eq!(info.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_file_format() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from(".")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert_eq!(
            Info::try_from_record_file_format(record, FileFormat::new(4, 2)),
            Ok(Info::new(
                Key::SamplesWithDataCount,
                Number::Unknown,
                Type::Integer,
                String::from("Number of samples with data"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_info() {
        let record = build_record();

        assert_eq!(
            Info::try_from(record),
            Ok(Info::new(
                Key::SamplesWithDataCount,
                Number::Count(1),
                Type::Integer,
                String::from("Number of samples with data"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_info_with_extra_fields() -> Result<(), &'static str> {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
                (String::from("Source"), String::from("dbsnp")),
                (String::from("Version"), String::from("138")),
                (String::from("IDX"), String::from("1")),
            ]),
        );

        assert_eq!(
            Info::try_from(record),
            Ok(Info {
                id: Key::SamplesWithDataCount,
                number: Number::Count(1),
                ty: Type::Integer,
                description: String::from("Number of samples with data"),
                idx: Some(1),
                fields: [
                    (
                        Tag::other("Source").ok_or("invalid tag")?,
                        String::from("dbsnp")
                    ),
                    (
                        Tag::other("Version").ok_or("invalid tag")?,
                        String::from("138")
                    ),
                ]
                .into_iter()
                .collect()
            })
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_record_key() {
        let record = Record::new(
            record::Key::FileFormat,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert_eq!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_record_value() {
        let record = Record::new(
            record::Key::Info,
            record::Value::String(String::from("VCFv4.3")),
        );

        assert_eq!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_info_with_a_missing_field() {
        let record = Record::new(record::Key::Info, record::Value::Struct(Vec::new()));

        assert_eq!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_id() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidId(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_number() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("NA")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidNumber(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_type() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("int")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidType(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_info_with_an_invalid_idx() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
                (String::from("IDX"), String::from("ndls")),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_info_with_a_number_mismatch() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from(".")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::NumberMismatch(
                Number::Unknown,
                Number::Count(1)
            ))
        ));
    }

    #[test]
    fn test_try_from_record_for_info_with_a_type_mismatch() {
        let record = Record::new(
            record::Key::Info,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("String")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data"),
                ),
            ]),
        );

        assert!(matches!(
            Info::try_from(record),
            Err(TryFromRecordError::TypeMismatch(
                Type::String,
                Type::Integer
            ))
        ));
    }
}
