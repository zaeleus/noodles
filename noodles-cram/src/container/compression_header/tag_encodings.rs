use std::collections::HashMap;

use super::{Encoding, encoding::codec::ByteArray};
use crate::container::block;

pub type TagEncodings = HashMap<block::ContentId, Encoding<ByteArray>>;
