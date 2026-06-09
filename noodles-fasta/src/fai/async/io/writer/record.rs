use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::fai::Record;

pub(super) async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::fai::io::writer::write_record;

    let mut buf = Vec::new();
    write_record(&mut buf, record)?;
    writer.write_all(&buf).await
}
