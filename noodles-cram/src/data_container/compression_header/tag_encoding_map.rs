mod builder;

pub use self::builder::Builder;

use std::collections::HashMap;

use super::{encoding::codec::ByteArray, Encoding};
use crate::container::block;

pub type TagEncodingMap = HashMap<block::ContentId, Encoding<ByteArray>>;
