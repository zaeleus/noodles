mod key;

use std::convert::TryFrom;

use self::key::Key;

#[derive(Clone, Debug)]
pub struct Filter {
    id: String,
    description: String,
}

impl Filter {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn description(&self) -> &str {
        &self.description
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
}

impl TryFrom<&[(String, String)]> for Filter {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut it = fields.iter();

        let id = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Id))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Id) => Ok(v.into()),
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
            (String::from("ID"), String::from("q10")),
            (
                String::from("Description"),
                String::from("Quality below 10"),
            ),
        ];

        let filter = Filter::try_from(&fields[..])?;

        assert_eq!(filter.id(), "q10");
        assert_eq!(filter.description(), "Quality below 10");

        Ok(())
    }
}
