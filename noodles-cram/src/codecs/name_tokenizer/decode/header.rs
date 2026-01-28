use std::io;

use crate::io::reader::num::{read_u8, read_u32_le};

pub(super) fn read_header(src: &mut &[u8]) -> io::Result<(usize, usize, bool)> {
    let ulen = read_u32_le(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_names = read_u32_le(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let use_arith = read_u8(src)? == 1;

    Ok((ulen, n_names, use_arith))
}
