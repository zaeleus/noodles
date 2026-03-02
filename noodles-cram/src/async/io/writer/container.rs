use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{
    file_definition::Version,
    io::writer::{Options, Record},
};

pub async fn write_container<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    options: &Options,
    header: &sam::Header,
    record_counter: u64,
    records: &mut [Record],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();

    crate::io::writer::container::write_container(
        &mut buf,
        reference_sequence_repository,
        options,
        header,
        record_counter,
        records,
    )?;

    writer.write_all(&buf).await?;

    Ok(())
}

pub async fn write_eof_container<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let eof = crate::io::writer::container::build_eof_container(version)?;
    writer.write_all(&eof).await
}
