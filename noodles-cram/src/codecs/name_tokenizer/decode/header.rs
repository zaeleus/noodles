use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_u32_le};

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<(usize, usize, bool)>
where
    R: Read,
{
    let ulen = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_names = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let use_arith = read_u8(reader)? == 1;

    Ok((ulen, n_names, use_arith))
}
