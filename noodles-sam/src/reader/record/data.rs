use std::{io, str};

use crate::record::Data;

pub(crate) fn parse_data(src: &[u8]) -> io::Result<Data> {
    if src.is_empty() {
        Ok(Data::default())
    } else {
        str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }
}
