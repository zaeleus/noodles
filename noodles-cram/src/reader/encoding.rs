use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{encoding, num::read_itf8, Encoding};

pub fn read_encoding<R>(reader: &mut R) -> io::Result<Encoding>
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
