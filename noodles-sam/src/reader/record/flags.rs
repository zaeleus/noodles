use std::io;

use crate::record::Flags;

pub(crate) fn parse_flags(src: &[u8]) -> io::Result<Flags> {
    lexical_core::parse::<u16>(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Flags::from)
}
