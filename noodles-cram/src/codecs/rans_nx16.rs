//! rANS 4x8 codec.

pub(crate) mod decode;
pub(crate) mod encode;
mod flags;

pub use self::flags::Flags;
pub(crate) use self::{decode::decode, encode::encode};
