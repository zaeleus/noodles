use std::{io, str};

use crate::Header;

pub(super) fn parse_reference_sequence_id(
    header: &Header,
    src: &[u8],
) -> io::Result<Option<usize>> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(None),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                header
                    .reference_sequences()
                    .get_index_of(s)
                    .map(Some)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid reference sequence name",
                        )
                    })
            }),
    }
}
