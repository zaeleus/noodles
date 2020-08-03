mod key;

pub use self::key::Key;

#[derive(Clone, Debug)]
pub struct Tag {
    key: Key,
    data: Vec<u8>,
}

impl Tag {
    pub fn new(key: Key, data: Vec<u8>) -> Self {
        Self { key, data }
    }
}
