mod builder;

pub use self::builder::Builder;

use std::{collections::HashMap, ops::Deref};

use super::Encoding;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagEncodingMap(HashMap<i32, Encoding>);

impl Deref for TagEncodingMap {
    type Target = HashMap<i32, Encoding>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<HashMap<i32, Encoding>> for TagEncodingMap {
    fn from(map: HashMap<i32, Encoding>) -> Self {
        Self(map)
    }
}
