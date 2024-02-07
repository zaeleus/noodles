use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::Header;

pub(super) async fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut serializer = crate::io::Writer::new(Vec::new());
    serializer.write_header(header)?;
    writer.write_all(serializer.get_ref()).await
}
