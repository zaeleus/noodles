mod builder;

pub use self::builder::Builder;

use std::{collections::HashMap, ops::Deref};

use super::{encoding::codec::ByteArray, Encoding};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagEncodingMap(HashMap<i32, Encoding<ByteArray>>);

impl Deref for TagEncodingMap {
    type Target = HashMap<i32, Encoding<ByteArray>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<HashMap<i32, Encoding<ByteArray>>> for TagEncodingMap {
    fn from(map: HashMap<i32, Encoding<ByteArray>>) -> Self {
        Self(map)
    }
}
