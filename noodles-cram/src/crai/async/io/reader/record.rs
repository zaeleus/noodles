use tokio::io::{self, AsyncBufRead};

use super::read_line;
use crate::crai::{Record, io::reader::record::parse_record};

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut String,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            *record =
                parse_record(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}
