use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use crate::lazy;

pub(super) async fn read_lazy_line<R>(
    reader: &mut R,
    buf: &mut String,
    line: &mut lazy::Line,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    buf.clear();

    if reader.read_line(buf).await? == 0 {
        return Ok(0);
    }

    let mut reader = buf.as_bytes();
    crate::io::reader::read_lazy_line(&mut reader, line)
}
