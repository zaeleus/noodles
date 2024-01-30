use std::io;

use byteorder::{LittleEndian, ReadBytesExt};

use crate::record::ChromosomeId;

pub(crate) fn read_chrom(src: &mut &[u8]) -> io::Result<ChromosomeId> {
    src.read_i32::<LittleEndian>().and_then(|n| {
        ChromosomeId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}
