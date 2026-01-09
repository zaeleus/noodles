#![doc(hidden)]

pub mod decoder;
pub mod encoder;

pub use self::encoder::{encode_record_buf, encode_with_prealloc, estimate_record_size};
pub(crate) use self::{decoder::decode, encoder::encode};
