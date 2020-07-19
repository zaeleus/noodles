use std::io::{self, Read};

use crate::{
    container::compression_header::TagEncodingMap, num::read_itf8, reader::encoding::read_encoding,
};

pub fn read_tag_encoding_map<R>(reader: &mut R) -> io::Result<TagEncodingMap>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut encodings = TagEncodingMap::with_capacity(map_len as usize);

    for _ in 0..map_len {
        let key = read_itf8(&mut buf_reader)?;
        let encoding = read_encoding(&mut buf_reader)?;
        encodings.insert(key, encoding);
    }

    Ok(encodings)
}
