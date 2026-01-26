use std::io;

use super::num::read_itf8_as;

pub(super) fn read_array<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8]> {
    let len = read_itf8_as(src)?;

    src.split_off(..len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

pub(super) fn read_map<'a>(src: &mut &'a [u8]) -> io::Result<(&'a [u8], usize)> {
    let mut buf = read_array(src)?;
    let len = read_itf8_as(&mut buf)?;
    Ok((buf, len))
}
