mod platform;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{platform::Platform, tag::Tag};

use super::record;

#[derive(Debug)]
pub struct ReadGroup {
    id: String,
    fields: HashMap<Tag, String>,
}

impl ReadGroup {
    pub fn new(id: String) -> Self {
        Self {
            id,
            ..Default::default()
        }
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }
}

impl Default for ReadGroup {
    fn default() -> Self {
        Self {
            id: String::new(),
            fields: HashMap::new(),
        }
    }
}

impl fmt::Display for ReadGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::ReadGroup)?;
        write!(f, "\t{}:{}", Tag::Id, self.id)?;

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
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

impl TryFrom<&[(String, String)]> for ReadGroup {
    type Error = ParseError;

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut read_group = ReadGroup::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

            if let Tag::Id = tag {
                read_group.id = value.into();
                has_id = true;
                continue;
            }

            read_group.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(ParseError::MissingRequiredTag(Tag::Id));
        }

        Ok(read_group)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let mut read_group = ReadGroup::new(String::from("rg0"));

        read_group
            .fields
            .insert(Tag::Program, String::from("noodles"));

        let actual = format!("{}", read_group);
        let expected = "@RG\tID:rg0\tPG:noodles";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_with_no_version() {
        let fields = [(String::from("DS"), String::from("noodles"))];
        assert!(ReadGroup::try_from(&fields[..]).is_err());
    }
}
