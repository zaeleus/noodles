mod tag;

use std::{collections::HashMap, convert::TryFrom};

pub use self::tag::Tag;

#[derive(Debug)]
pub struct Program(HashMap<Tag, String>);

impl Program {
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.0.get(tag)
    }
}

impl Default for Program {
    fn default() -> Self {
        Self(HashMap::new())
    }
}

impl TryFrom<&[(String, String)]> for Program {
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
