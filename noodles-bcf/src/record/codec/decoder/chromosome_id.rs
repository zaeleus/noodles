use std::io;

use byteorder::{LittleEndian, ReadBytesExt};

pub(crate) fn read_chrom(src: &mut &[u8]) -> io::Result<usize> {
    src.read_i32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}
