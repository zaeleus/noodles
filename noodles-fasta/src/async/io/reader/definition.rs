use tokio::io::{self, AsyncBufRead};

use super::read_line;
use crate::record::Definition;

pub(super) async fn read_definition<R>(
    reader: &mut R,
    buf: &mut String,
    definition: &mut Definition,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    buf.clear();

    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            *definition = buf
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}
