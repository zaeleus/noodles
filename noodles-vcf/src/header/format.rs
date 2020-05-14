mod key;
mod ty;

pub use self::ty::Type;

use std::convert::TryFrom;

use crate::record::format;

use super::{number, Number};

use self::key::Key;

#[derive(Clone, Debug)]
pub struct Format {
    id: format::Key,
    number: Number,
    ty: Type,
    description: String,
}

impl Format {
    pub fn id(&self) -> &format::Key {
        &self.id
    }

    pub fn number(&self) -> &Number {
        &self.number
    }

    pub fn ty(&self) -> &Type {
        &self.ty
    }

    pub fn description(&self) -> &str {
        &self.description
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
    InvalidId(format::key::ParseError),
    InvalidNumber(number::ParseError),
    InvalidType(ty::ParseError),
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

    #[test]
    fn test_try_from_fields_for_format() -> Result<(), ParseError> {
        let fields = vec![
            (String::from("ID"), String::from("GT")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (String::from("Description"), String::from("Genotype")),
        ];

        let format = Format::try_from(&fields[..])?;

        assert_eq!(format.id(), &format::Key::Genotype);
        assert_eq!(*format.number(), Number::Count(1));
        assert_eq!(*format.ty(), Type::Integer);
        assert_eq!(format.description(), "Genotype");

        Ok(())
    }
}
