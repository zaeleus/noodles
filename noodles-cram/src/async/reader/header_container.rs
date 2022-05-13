mod header;

use bytes::BytesMut;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::header::read_header;

pub async fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    use crate::reader::header_container::read_raw_sam_header_from_block;

    let len = read_header(reader).await?;

    buf.resize(len, 0);
    reader.read_exact(buf).await?;
    let mut buf = buf.split().freeze();

    read_raw_sam_header_from_block(&mut buf)
}
