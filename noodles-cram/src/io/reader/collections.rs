use std::io;

use crate::file_definition::Version;

use super::num::read_unsigned_int_as;

pub(super) fn read_array<'a>(src: &mut &'a [u8], version: Version) -> io::Result<&'a [u8]> {
    let len = read_unsigned_int_as(src, version)?;

    src.split_off(..len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

pub(super) fn read_map<'a>(src: &mut &'a [u8], version: Version) -> io::Result<(&'a [u8], usize)> {
    let mut buf = read_array(src, version)?;
    let len = read_unsigned_int_as(&mut buf, version)?;
    Ok((buf, len))
}
