use std::io;

use bytes::{Buf, Bytes};

use super::get_encoding_for_byte_array_codec;
use crate::{container::compression_header::TagEncodings, io::reader::num::get_itf8};

pub fn get_tag_encodings(src: &mut Bytes) -> io::Result<TagEncodings> {
    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let map_len = get_itf8(&mut buf).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut encodings = TagEncodings::with_capacity(map_len);

    for _ in 0..map_len {
        let key = get_itf8(&mut buf)?;
        let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
        encodings.insert(key, encoding);
    }

    Ok(encodings)
}
