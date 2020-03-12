mod group_order;
mod sort_order;
mod subsort_order;
mod tag;

use std::{collections::HashMap, convert::TryFrom};

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, tag::Tag,
};

#[derive(Debug)]
pub struct Header {
    version: String,
    fields: HashMap<Tag, String>,
}

impl Header {
    pub fn version(&self) -> &str {
        &self.version
    }

    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }
}

impl Default for Header {
    fn default() -> Self {
        Header {
            version: String::new(),
            fields: HashMap::new(),
        }
    }
}

impl TryFrom<&[(String, String)]> for Header {
    type Error = ();

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut header = Header::default();

        let mut has_version = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse()?;

            if let Tag::Version = tag {
                header.version = value.into();
                has_version = true;
            }

            header.fields.insert(tag, value.into());
        }

        if !has_version {
            return Err(());
        }

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_with_no_version() {
        let fields = [(String::from("SO"), String::from("coordinate"))];
        assert!(Header::try_from(&fields[..]).is_err());
    }
}
