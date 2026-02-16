use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::FileDefinition;

pub(super) async fn write_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    file_definition: &FileDefinition,
    header: &sam::Header,
    reference_required: bool,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_file_definition(writer, file_definition).await?;
    write_file_header(
        writer,
        reference_sequence_repository,
        header,
        file_definition.version(),
        reference_required,
    )
    .await?;
    Ok(())
}

pub(super) async fn write_file_definition<W>(
    writer: &mut W,
    file_definition: &FileDefinition,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::header::write_file_definition(&mut buf, file_definition)?;
    writer.write_all(&buf).await?;
    Ok(())
}

pub(super) async fn write_file_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    version: crate::file_definition::Version,
    reference_required: bool,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mut buf = Vec::new();
    crate::io::writer::header::write_file_header(
        &mut buf,
        reference_sequence_repository,
        header,
        version,
        reference_required,
    )?;
    writer.write_all(&buf).await?;
    Ok(())
}
