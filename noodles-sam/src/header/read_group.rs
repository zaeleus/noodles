mod platform;
mod tag;

use std::{collections::HashMap, convert::TryFrom};

pub use self::{platform::Platform, tag::Tag};

#[derive(Debug)]
pub struct ReadGroup(HashMap<Tag, String>);

impl ReadGroup {
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.0.get(tag)
    }
}

impl Default for ReadGroup {
    fn default() -> Self {
        Self(HashMap::new())
    }
}

impl TryFrom<&[(String, String)]> for ReadGroup {
    type Error = ();

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut fields = HashMap::new();

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse()?;
            fields.insert(tag, value.into());
        }

        Ok(Self(fields))
    }
}
