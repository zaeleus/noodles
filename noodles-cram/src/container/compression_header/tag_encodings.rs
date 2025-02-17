use std::collections::HashMap;

use super::{encoding::codec::ByteArray, Encoding};
use crate::container::block;

pub type TagEncodings = HashMap<block::ContentId, Encoding<ByteArray>>;
