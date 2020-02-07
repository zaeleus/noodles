use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{encoding, num::read_itf8, Encoding, TagEncodings};

pub fn read_tag_encodings<R>(reader: &mut R) -> io::Result<TagEncodings>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut encodings = TagEncodings::with_capacity(map_len as usize);

    for _ in 0..map_len {
        let key = read_itf8(&mut buf_reader)?;
        let encoding = read_encoding(&mut buf_reader)?;
        encodings.insert(key, encoding);
    }

    Ok(encodings)
}

fn read_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let kind = read_itf8(reader)
        .map(|codec_id| encoding::Kind::try_from(codec_id).expect("invalid codec id"))?;

    let args_len = read_itf8(reader)?;
    let mut args_buf = vec![0; args_len as usize];
    reader.read_exact(&mut args_buf)?;

    Ok(Encoding::new(kind, args_buf))
}
