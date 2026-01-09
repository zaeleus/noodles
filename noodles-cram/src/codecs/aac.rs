//! Adaptive arithmetic coder.

mod decode;
mod encode;
mod flags;
mod model;
mod range_coder;
mod rle;

pub use self::flags::Flags;
pub(crate) use self::{decode::decode, encode::encode, model::Model, range_coder::RangeCoder};
