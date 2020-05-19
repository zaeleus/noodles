use noodles_sam::record::data::field::Tag;

use super::Value;

#[derive(Debug)]
pub struct Field {
    tag: Tag,
    value: Value,
}

impl Field {
    pub fn new(tag: Tag, value: Value) -> Self {
        Self { tag, value }
    }

    pub fn tag(&self) -> &Tag {
        &self.tag
    }

    pub fn value(&self) -> &Value {
        &self.value
    }
}
