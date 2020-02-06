use std::io::{self, Read};

use crate::num::read_itf8;

pub fn read_tag_encodings<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;

    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    for _ in 0..map_len {
        let _key = read_itf8(&mut buf_reader)?;
        read_encoding(&mut buf_reader)?
    }

    Ok(())
}

fn read_encoding<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let _codec_id = read_itf8(reader)?;
    let len = read_itf8(reader)?;
    let mut buf = vec![0; len as usize];
    reader.read_exact(&mut buf)
}
