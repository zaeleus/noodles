use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::Position;

pub fn read_pos(src: &mut &[u8]) -> io::Result<Position> {
    src.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n + 1)
            .map(Position::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}
