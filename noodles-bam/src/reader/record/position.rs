use std::{io, mem};

use bytes::Buf;
use noodles_core::Position;

pub(crate) fn get_position<B>(src: &mut B) -> io::Result<Option<Position>>
where
    B: Buf,
{
    const MISSING: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_i32_le() {
        MISSING => Ok(None),
        n => usize::try_from(n + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Position::new),
    }
}
