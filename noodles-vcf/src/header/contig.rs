mod key;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

use super::record;

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

impl fmt::Display for Contig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Kind::Contig.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;

        for (key, value) in &self.fields {
            write!(f, r#",{}="{}""#, key, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

#[derive(Debug)]
pub enum ParseError {
    InvalidKey(key::ParseError),
    MissingField(Key),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid contig header: ")?;

        match self {
            ParseError::InvalidKey(e) => write!(f, "{}", e),
            ParseError::MissingField(key) => write!(f, "missing {} field", key),
        }
    }
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

    fn build_fields() -> Vec<(String, String)> {
        vec![
            (String::from("ID"), String::from("sq0")),
            (String::from("length"), String::from("13")),
        ]
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let fields = build_fields();
        let contig = Contig::try_from(&fields[..])?;

        let expected = r#"##contig=<ID=sq0,length="13">"#;

        assert_eq!(contig.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_contig() -> Result<(), ParseError> {
        let fields = build_fields();
        let contig = Contig::try_from(&fields[..])?;

        assert_eq!(contig.id(), "sq0");
        assert_eq!(contig.get("length"), Some("13"));

        Ok(())
    }
}
