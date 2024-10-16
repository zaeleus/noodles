mod definition;
mod sequence;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{definition::write_definition, sequence::write_sequence};
use crate::Record;

const LINE_BASES: usize = 80;

pub(super) async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_definition(writer, record.definition()).await?;
    write_newline(writer).await?;

    write_sequence(writer, record.sequence(), LINE_BASES).await?;

    Ok(())
}

async fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    writer.write_all(&[LINE_FEED]).await
}
