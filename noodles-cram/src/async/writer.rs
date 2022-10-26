use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{file_definition::Version, FileDefinition, MAGIC_NUMBER};

/// An async CRAM writer.
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
