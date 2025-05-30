use tokio::io::{self, AsyncBufRead};

use super::read_line;
use crate::{Header, alignment::RecordBuf};

pub(super) async fn read_record_buf<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    header: &Header,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use crate::io::reader::record_buf::parse_record_buf;

    buf.clear();

    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            parse_record_buf(buf, header, record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}
