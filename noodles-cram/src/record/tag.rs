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
}
