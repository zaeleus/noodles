mod key;
mod ty;

pub use self::ty::Type;

use std::{convert::TryFrom, error, fmt};

use crate::record::genotype;

use super::{number, record, Number};

use self::key::Key;

#[derive(Clone, Debug)]
pub struct Format {
    id: genotype::field::Key,
    number: Number,
    ty: Type,
    description: String,
}

impl Format {
    pub fn new(id: genotype::field::Key, number: Number, ty: Type, description: String) -> Self {
        Self {
            id,
            number,
            ty,
            description,
        }
    }

    pub fn id(&self) -> &genotype::field::Key {
        &self.id
    }

    pub fn number(&self) -> Number {
        self.number
    }

    pub fn ty(&self) -> Type {
        self.ty
    }

    pub fn description(&self) -> &str {
        &self.description
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Kind::Format.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, ",{}={}", Key::Number, self.number)?;
        write!(f, ",{}={}", Key::Type, self.ty)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

        f.write_str(">")?;

        Ok(())
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
    InvalidId(genotype::field::key::ParseError),
    InvalidNumber(number::ParseError),
    InvalidType(ty::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid format header: ")?;

        match self {
            ParseError::MissingField(key) => write!(f, "missing {} field", key),
            ParseError::InvalidId(e) => write!(f, "{}", e),
            ParseError::InvalidNumber(e) => write!(f, "{}", e),
            ParseError::InvalidType(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<&[(String, String)]> for Format {
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
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_fields() -> Vec<(String, String)> {
        vec![
            (String::from("ID"), String::from("GT")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (String::from("Description"), String::from("Genotype")),
        ]
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let fields = build_fields();
        let format = Format::try_from(&fields[..])?;

        let expected = r#"##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">"#;

        assert_eq!(format.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_format() -> Result<(), ParseError> {
        let fields = build_fields();

        let format = Format::try_from(&fields[..])?;

        assert_eq!(format.id(), &genotype::field::Key::Genotype);
        assert_eq!(format.number(), Number::Count(1));
        assert_eq!(format.ty(), Type::Integer);
        assert_eq!(format.description(), "Genotype");

        Ok(())
    }
}
