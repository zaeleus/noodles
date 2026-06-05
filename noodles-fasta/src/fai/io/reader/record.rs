use std::io::{self, BufRead};

use super::read_line;
use crate::fai::Record;

pub(super) fn read_record<R>(
    reader: &mut R,
    buf: &mut String,
    record: &mut Record,
) -> io::Result<usize>
where
    R: BufRead,
{
    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            *record = buf
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}
