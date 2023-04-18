use std::io;

use crate::record::ReadName;

pub(crate) fn parse_read_name(src: &[u8]) -> io::Result<ReadName> {
    ReadName::try_new(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
