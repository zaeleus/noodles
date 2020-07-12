mod key;
mod ty;

pub use self::ty::Type;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

use crate::record::info;

use super::{number, record, Number};

use self::key::Key;

/// A VCF header information record (`INFO`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Info {
    id: info::field::Key,
    number: Number,
    ty: Type,
    description: String,
    fields: HashMap<String, String>,
}

impl Info {
    /// Creates a VCF header information record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
    /// );
    /// ```
    pub fn new(id: info::field::Key, number: Number, ty: Type, description: String) -> Self {
        Self {
            id,
            number,
            ty,
            description,
            fields: HashMap::new(),
        }
    }

    /// Returns the information field key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
    /// );
    ///
    /// assert_eq!(info.id(), &Key::SamplesWithDataCount);
    /// ```
    pub fn id(&self) -> &info::field::Key {
        &self.id
    }

    /// Returns the cardinality of the information field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
    /// );
    ///
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
    /// use noodles_vcf::{
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
    /// );
    ///
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
    /// use noodles_vcf::{
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let info = Info::new(
    ///     Key::SamplesWithDataCount,
    ///     Number::Count(1),
    ///     Type::Integer,
    ///     String::from("Number of samples with data"),
    /// );
    ///
    /// assert_eq!(info.description(), "Number of samples with data");
    /// ```
    pub fn description(&self) -> &str {
        &self.description
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID`, `Number`, `Type`, and `Description`.
    pub fn fields(&self) -> &HashMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Key::Info.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, ",{}={}", Key::Number, self.number)?;
        write!(f, ",{}={}", Key::Type, self.ty)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

        for (key, value) in &self.fields {
            write!(f, r#",{}="{}""#, key, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a raw VCF header information record fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// A required field is missing.
    MissingField(Key),
    /// The ID is invalid.
    InvalidId(info::field::key::ParseError),
    /// The number is invalid.
    InvalidNumber(number::ParseError),
    /// The type is invalid.
    InvalidType(ty::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid info header: ")?;

        match self {
            ParseError::MissingField(key) => write!(f, "missing {} field", key),
            ParseError::InvalidId(e) => write!(f, "{}", e),
            ParseError::InvalidNumber(e) => write!(f, "{}", e),
            ParseError::InvalidType(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<&[(String, String)]> for Info {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut it = fields.iter();

        let id = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Id))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Id) => v.parse().map_err(ParseError::InvalidId),
                _ => Err(ParseError::MissingField(Key::Id)),
            })?;

        let number = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Number))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Number) => v.parse().map_err(ParseError::InvalidNumber),
                _ => Err(ParseError::MissingField(Key::Id)),
            })?;

        let ty = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Type))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Type) => v.parse().map_err(ParseError::InvalidType),
                _ => Err(ParseError::MissingField(Key::Type)),
            })?;

        let description = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Description))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Description) => Ok(v.into()),
                _ => Err(ParseError::MissingField(Key::Description)),
            })?;

        Ok(Self {
            id,
            number,
            ty,
            description,
            fields: it.cloned().collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_fields() -> Vec<(String, String)> {
        vec![
            (String::from("ID"), String::from("NS")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (
                String::from("Description"),
                String::from("Number of samples with data"),
            ),
        ]
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let fields = build_fields();
        let info = Info::try_from(&fields[..])?;

        let expected =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;

        assert_eq!(info.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_info() -> Result<(), ParseError> {
        let fields = build_fields();
        let info = Info::try_from(&fields[..])?;

        assert_eq!(info.id(), &info::field::Key::SamplesWithDataCount);
        assert_eq!(info.number(), Number::Count(1));
        assert_eq!(info.ty(), Type::Integer);
        assert_eq!(info.description(), "Number of samples with data");

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_info_with_extra_fields() -> Result<(), ParseError> {
        let mut fields = build_fields();
        fields.push((String::from("Source"), String::from("dbsnp")));
        fields.push((String::from("Version"), String::from("138")));

        let info = Info::try_from(&fields[..])?;

        let fields = info.fields();
        assert_eq!(fields.len(), 2);
        assert_eq!(fields["Source"], "dbsnp");
        assert_eq!(fields["Version"], "138");

        Ok(())
    }
}
