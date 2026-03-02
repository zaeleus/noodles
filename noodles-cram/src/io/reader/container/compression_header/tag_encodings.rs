use std::io;

use super::read_byte_array_encoding;
use crate::{
    container::compression_header::TagEncodings,
    file_definition::Version,
    io::reader::{collections::read_map, num::read_unsigned_int},
};

pub fn read_tag_encodings(src: &mut &[u8], version: Version) -> io::Result<TagEncodings> {
    let (mut buf, len) = read_map(src, version)?;
    read_tag_encodings_inner(&mut buf, len, version)
}

fn read_tag_encodings_inner(
    src: &mut &[u8],
    len: usize,
    version: Version,
) -> io::Result<TagEncodings> {
    (0..len)
        .map(|_| {
            let block_content_id = read_unsigned_int(src, version)?;
            let encoding = read_byte_array_encoding(src, version)?;
            Ok((block_content_id, encoding))
        })
        .collect()
}
