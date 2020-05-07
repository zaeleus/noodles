mod key;

use std::convert::TryFrom;

use self::key::Key;

#[derive(Clone, Debug)]
pub struct Contig {
    id: String,
}

impl Contig {
    pub fn id(&self) -> &str {
        &self.id
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
}

impl TryFrom<&[(String, String)]> for Contig {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut it = fields.iter();

        let id = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Id))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Id) => Ok(v.into()),
                _ => return Err(ParseError::MissingField(Key::Id)),
            })?;

        Ok(Self { id })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_fields_for_contig() -> Result<(), ParseError> {
        let fields = vec![(String::from("ID"), String::from("sq0"))];
        let contig = Contig::try_from(&fields[..])?;

        assert_eq!(contig.id(), "sq0");

        Ok(())
    }
}
