use crate::{
    container::block,
    data_container::compression_header::{
        encoding::codec::{Byte, Integer},
        Encoding,
    },
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ByteArray {
    // len_encoding, value_encoding
    ByteArrayLen(Encoding<Integer>, Encoding<Byte>),
    // stop_byte, block_content_id
    ByteArrayStop(u8, block::ContentId),
}
