use std::{io, str};

use crate::record::Sequence;

pub(crate) fn parse_sequence(src: &[u8]) -> io::Result<Sequence> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(Sequence::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}
