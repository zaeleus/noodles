use std::io;

use crate::{Header, Record};

pub(super) fn parse_record(src: &str, header: &Header, record: &mut Record) -> io::Result<()> {
    *record = Record::try_from_str(src, header)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(())
}
