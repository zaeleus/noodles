mod tag;

use std::{collections::HashMap, convert::TryFrom};

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

impl TryFrom<&[(String, String)]> for Program {
    type Error = ();

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut program = Program::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse()?;

            if let Tag::Id = tag {
                program.id = value.into();
                has_id = true;
            }

            program.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(());
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
