use crate::{r#async::io::reader::read_line, fai::Record};
use tokio::io::{self, AsyncBufRead};

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
            *record = buf
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}
