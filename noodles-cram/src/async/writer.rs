mod container;
mod header_container;

use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{file_definition::Version, FileDefinition, MAGIC_NUMBER};

/// An async CRAM writer.
///
/// A call to [`Self::shutdown`] must be made before the writer is dropped.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::AsyncWriter::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::AsyncWriter::new(io::sink());
    /// let inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Attempts to shutdown the output stream by writing any pending containers and a final EOF
    /// container.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let mut writer = cram::AsyncWriter::new(io::sink());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        container::write_eof_container(&mut self.inner).await
    }

    /// Writes a CRAM file definition.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_cram as cram;
    ///
    /// let mut writer = cram::AsyncWriter::new(Vec::new());
    /// writer.write_file_definition().await?;
    ///
    /// assert_eq!(writer.get_ref(), &[
    ///     // magic number (CRAM)
    ///     0x43, 0x52, 0x41, 0x4d,
    ///     // format (major, minor)
    ///     0x03, 0x00,
    ///     // file ID
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    /// ]);
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_file_definition(&mut self) -> io::Result<()> {
        let file_definition = FileDefinition::default();
        write_file_definition(&mut self.inner, &file_definition).await
    }

    /// Writes a CRAM file header container.
    ///
    /// The position of the stream is expected to be directly after the file definition.
    ///
    /// Entries in the reference sequence dictionary that are missing MD5 checksums (`M5`) will
    /// automatically be calculated and added to the written record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = cram::AsyncWriter::new(io::sink());
    /// writer.write_file_definition().await?;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// writer.write_file_header(&repository, &header).await?;
    ///
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_file_header(
        &mut self,
        reference_sequence_repository: &fasta::Repository,
        header: &sam::Header,
    ) -> io::Result<()> {
        use self::header_container::write_header_container;
        use crate::writer::add_missing_reference_sequence_checksums;

        let mut header = header.clone();

        add_missing_reference_sequence_checksums(
            reference_sequence_repository,
            header.reference_sequences_mut(),
        )?;

        write_header_container(&mut self.inner, &header).await
    }
}

async fn write_file_definition<W>(
    writer: &mut W,
    file_definition: &FileDefinition,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(MAGIC_NUMBER).await?;
    write_format(writer, file_definition.version()).await?;
    writer.write_all(file_definition.file_id()).await?;
    Ok(())
}

async fn write_format<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let format = [version.major(), version.minor()];
    writer.write_all(&format).await
}
