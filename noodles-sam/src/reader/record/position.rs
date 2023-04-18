use std::io;

use noodles_core::Position;

pub(crate) fn parse_alignment_start(src: &[u8]) -> io::Result<Option<Position>> {
    lexical_core::parse(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Position::new)
}
