mod platform;
mod tag;

use std::{collections::HashMap, convert::TryFrom};

pub use self::{platform::Platform, tag::Tag};

#[derive(Debug)]
pub struct ReadGroup {
    id: String,
    fields: HashMap<Tag, String>,
}

impl ReadGroup {
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

impl TryFrom<&[(String, String)]> for ReadGroup {
    type Error = ();

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut read_group = ReadGroup::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse()?;

            if let Tag::Id = tag {
                read_group.id = value.into();
                has_id = true;
            }

            read_group.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(());
        }

        Ok(read_group)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_with_no_version() {
        let fields = [(String::from("DS"), String::from("noodles"))];
        assert!(ReadGroup::try_from(&fields[..]).is_err());
    }
}
