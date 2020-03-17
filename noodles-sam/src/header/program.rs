mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::tag::Tag;

#[derive(Debug)]
pub struct Program {
    id: String,
    fields: HashMap<Tag, String>,
}

impl Program {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }
}

impl Default for Program {
    fn default() -> Self {
        Self {
            id: String::new(),
            fields: HashMap::new(),
        }
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingRequiredTag(Tag),
    InvalidTag(tag::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<&[(String, String)]> for Program {
    type Error = ParseError;

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut program = Program::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

            if let Tag::Id = tag {
                program.id = value.into();
                has_id = true;
            }

            program.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(ParseError::MissingRequiredTag(Tag::Id));
        }

        Ok(program)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_with_no_id() {
        let fields = [(String::from("DS"), String::from("noodles"))];
        assert!(Program::try_from(&fields[..]).is_err());
    }
}
