//! rANS Nx16 codec.

pub(crate) mod decode;
pub(crate) mod encode;
mod flags;

pub use self::flags::Flags;
pub(crate) use self::{decode::decode, encode::encode};

const ALPHABET_SIZE: usize = 256;
