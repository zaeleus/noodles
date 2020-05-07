mod key;

use std::{collections::HashMap, convert::TryFrom};

use self::key::Key;

#[derive(Clone, Debug)]
pub struct Contig {
    id: String,
    fields: HashMap<String, String>,
}

impl Contig {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        self.fields.get(key).map(|s| &**s)
    }
}

#[derive(Debug)]
pub enum ParseError {
    InvalidKey(key::ParseError),
    MissingField(Key),
}

impl TryFrom<&[(String, String)]> for Contig {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut contig = Self {
            id: String::from("unknown"),
            fields: HashMap::new(),
        };

        let mut has_id = false;

        for (raw_key, value) in fields {
            let key = raw_key.parse().map_err(ParseError::InvalidKey)?;

            match key {
                Key::Id => {
                    contig.id = value.into();
                    has_id = true;
                }
                Key::Other(k) => {
                    contig.fields.insert(k, value.into());
                }
            }
        }

        if !has_id {
            return Err(ParseError::MissingField(Key::Id));
        }

        Ok(contig)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_fields_for_contig() -> Result<(), ParseError> {
        let fields = vec![
            (String::from("ID"), String::from("sq0")),
            (String::from("length"), String::from("13")),
        ];

        let contig = Contig::try_from(&fields[..])?;

        assert_eq!(contig.id(), "sq0");
        assert_eq!(contig.get("length"), Some("13"));

        Ok(())
    }
}
