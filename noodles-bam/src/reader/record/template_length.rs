use std::{io, mem};

use bytes::Buf;

pub(super) fn get_template_length<B>(src: &mut B) -> io::Result<i32>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(src.get_i32_le())
}
