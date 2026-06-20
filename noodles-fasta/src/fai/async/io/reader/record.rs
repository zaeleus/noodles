use crate::{r#async::io::reader::read_line, fai::Record, fai::io::reader::record::parse_record};
use tokio::io::{self, AsyncBufRead};

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            let s =
                str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            *record = parse_record(s)?;

            Ok(n)
        }
    }
}
