//! BAM record codec for encoding and decoding.

pub mod decoder;
pub mod encoder;

pub use self::decoder::{DecodeError, decode};
pub use self::encoder::encode;
