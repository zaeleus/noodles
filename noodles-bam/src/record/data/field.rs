mod value;

pub use self::value::Value;

use noodles_sam::record::data::field::Tag;

#[derive(Clone, Debug, PartialEq)]
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
