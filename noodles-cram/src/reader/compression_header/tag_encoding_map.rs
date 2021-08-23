use std::{
    collections::HashMap,
    io::{self, Read},
};

use crate::{
    data_container::compression_header::TagEncodingMap,
    reader::{encoding::read_encoding, num::read_itf8},
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

    let mut map = HashMap::with_capacity(map_len as usize);

    for _ in 0..map_len {
        let key = read_itf8(&mut buf_reader)?;
        let encoding = read_encoding(&mut buf_reader)?;
        map.insert(key, encoding);
    }

    Ok(TagEncodingMap::from(map))
}
