mod header;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::header::read_header;
use crate::io::reader::Container;

pub async fn read_container<R>(reader: &mut R, container: &mut Container) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    match read_header(reader, &mut container.header).await? {
        0 => Ok(0),
        len => {
            container.src.resize(len, 0);
            reader.read_exact(&mut container.src).await?;
            Ok(len)
        }
    }
}
