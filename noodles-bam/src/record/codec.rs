pub mod decoder;
pub mod encoder;

pub(crate) use self::{decoder::decode_record as decode, encoder::encode_record as encode};
