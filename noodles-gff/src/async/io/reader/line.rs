use tokio::io::{self, AsyncBufRead};

use crate::Line;

pub(super) async fn read_line<R>(reader: &mut R, line: &mut Line) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    let buf = &mut line.0;
    buf.clear();
    super::read_line(reader, buf).await
}
