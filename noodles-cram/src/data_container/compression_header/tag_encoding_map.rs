mod builder;

pub use self::builder::Builder;

use std::{collections::HashMap, ops::Deref};

use super::{encoding::codec::ByteArray, Encoding};
use crate::container::block;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagEncodingMap(HashMap<block::ContentId, Encoding<ByteArray>>);

impl Deref for TagEncodingMap {
    type Target = HashMap<block::ContentId, Encoding<ByteArray>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<HashMap<block::ContentId, Encoding<ByteArray>>> for TagEncodingMap {
    fn from(map: HashMap<block::ContentId, Encoding<ByteArray>>) -> Self {
        Self(map)
    }
}
