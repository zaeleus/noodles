mod key;

pub use self::key::Key;

use noodles_bam::record::data::field::Value;

#[derive(Clone, Debug)]
pub struct Tag {
    key: Key,
    value: Value,
}

impl Tag {
    pub fn new(key: Key, value: Value) -> Self {
        Self { key, value }
    }

    pub fn key(&self) -> Key {
        self.key
    }

    pub fn value(&self) -> &Value {
        &self.value
    }
}
