use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::Flags;

pub(crate) fn get_flags<B>(src: &mut B) -> io::Result<Flags>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Flags::from(src.get_u16_le()))
}
