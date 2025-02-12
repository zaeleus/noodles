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

    let len = get_itf8(&mut buf).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    get_tag_encodings_inner(&mut buf, len)
}

fn get_tag_encodings_inner(src: &mut Bytes, len: usize) -> io::Result<TagEncodings> {
    let mut encodings = TagEncodings::with_capacity(len);

    for _ in 0..len {
        let key = get_itf8(src)?;
        let encoding = get_encoding_for_byte_array_codec(src)?;
        encodings.insert(key, encoding);
    }

    Ok(encodings)
}
