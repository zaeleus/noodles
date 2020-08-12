mod builder;

pub use self::builder::Builder;

use std::{collections::HashMap, ops::Deref};

use crate::num::Itf8;

use super::Encoding;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagEncodingMap(HashMap<Itf8, Encoding>);

impl TagEncodingMap {
    pub fn builder() -> Builder {
        Builder::default()
    }
}

impl Deref for TagEncodingMap {
    type Target = HashMap<Itf8, Encoding>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<HashMap<Itf8, Encoding>> for TagEncodingMap {
    fn from(map: HashMap<Itf8, Encoding>) -> Self {
        Self(map)
    }
}
