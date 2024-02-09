use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_core::Position;

pub fn read_pos(src: &mut &[u8]) -> io::Result<Option<Position>> {
    const TELOMERE_START: i32 = -1;

    match src.read_i32::<LittleEndian>()? {
        TELOMERE_START => Ok(None),
        n => usize::try_from(n + 1)
            .and_then(Position::try_from)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}
