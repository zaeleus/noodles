use std::{io, num::NonZeroUsize};

use bytes::{Buf, BytesMut};

use crate::record::ReadName;

pub fn get_read_name(
    src: &mut BytesMut,
    read_name: &mut Option<ReadName>,
    l_read_name: NonZeroUsize,
) -> io::Result<()> {
    const NUL: u8 = 0x00;
    const MISSING: [u8; 2] = [b'*', NUL];

    let len = usize::from(l_read_name);

    if src.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    *read_name = if src.take(len).chunk() == MISSING {
        src.advance(MISSING.len());
        None
    } else {
        let mut name = ReadName::default();
        name.buf = src.split_to(len - 1);
        src.advance(1);
        Some(name)
    };

    Ok(())
}
