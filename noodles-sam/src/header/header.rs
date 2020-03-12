mod group_order;
mod sort_order;
mod subsort_order;
mod tag;

use std::{collections::HashMap, convert::TryFrom};

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, tag::Tag,
};

#[derive(Debug)]
pub struct Header(HashMap<Tag, String>);

impl Header {
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.0.get(tag)
    }
}

impl Default for Header {
    fn default() -> Self {
        Self(HashMap::new())
    }
}

impl TryFrom<&[(String, String)]> for Header {
    type Error = ();

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut fields = HashMap::new();

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse()?;
            fields.insert(tag, value.into());
        }

        Ok(Header(fields))
    }
}
