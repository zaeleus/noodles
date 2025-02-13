use std::io;

use super::read_byte_array_encoding;
use crate::{
    container::compression_header::TagEncodings,
    io::reader::{collections::read_map, num::read_itf8},
};

pub fn read_tag_encodings(src: &mut &[u8]) -> io::Result<TagEncodings> {
    let (mut buf, len) = read_map(src)?;
    read_tag_encodings_inner(&mut buf, len)
}

fn read_tag_encodings_inner(src: &mut &[u8], len: usize) -> io::Result<TagEncodings> {
    (0..len)
        .map(|_| {
            let block_content_id = read_itf8(src)?;
            let encoding = read_byte_array_encoding(src)?;
            Ok((block_content_id, encoding))
        })
        .collect()
}
