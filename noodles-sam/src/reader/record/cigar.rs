use std::{io, str};

use crate::record::Cigar;

pub(crate) fn parse_cigar(src: &[u8]) -> io::Result<Cigar> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(Cigar::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}
