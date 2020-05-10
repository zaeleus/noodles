mod id;
mod key;

pub use self::id::Id;

use std::convert::TryFrom;

use self::key::Key;

#[derive(Clone, Debug)]
pub struct AlternativeAllele {
    id: Id,
    description: String,
}

impl AlternativeAllele {
    pub fn id(&self) -> Id {
        self.id
    }

    pub fn description(&self) -> &str {
        &self.description
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
    InvalidId(id::ParseError),
}

impl TryFrom<&[(String, String)]> for AlternativeAllele {
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

        let description = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Description))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Description) => Ok(v.into()),
                _ => Err(ParseError::MissingField(Key::Description)),
            })?;

        Ok(Self { id, description })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_fields_for_filter() -> Result<(), ParseError> {
        let fields = vec![
            (String::from("ID"), String::from("DEL")),
            (String::from("Description"), String::from("Deletion")),
        ];

        let filter = AlternativeAllele::try_from(&fields[..])?;

        assert_eq!(filter.id(), Id::Deletion);
        assert_eq!(filter.description(), "Deletion");

        Ok(())
    }
}
