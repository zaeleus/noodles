//! Async CRAM writer.

mod builder;
mod container;
mod data_container;
mod header_container;

use std::mem;

use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub use self::builder::Builder;
use crate::{
    file_definition::Version, io::writer::Options, DataContainer, FileDefinition, Record,
    MAGIC_NUMBER,
};

/// An async CRAM writer.
///
/// A call to [`Self::shutdown`] must be made before the writer is dropped.
pub struct Writer<W> {
    inner: W,
    reference_sequence_repository: fasta::Repository,
    options: Options,
    data_container_builder: crate::data_container::Builder,
    record_counter: u64,
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
    /// let writer = cram::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_from_writer(inner)
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::r#async::io::Writer::new(io::sink());
    /// let inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
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
    /// use noodles_sam as sam;
    /// use tokio::io;
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// let header = sam::Header::default();
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::container::write_eof_container;
        self.flush(header).await?;
        write_eof_container(&mut self.inner).await
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
    /// let mut writer = cram::r#async::io::Writer::new(Vec::new());
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
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// writer.write_file_definition().await?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header).await?;
    ///
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::header_container::write_header_container;
        use crate::io::writer::add_missing_reference_sequence_checksums;

        let mut header = header.clone();

        add_missing_reference_sequence_checksums(
            &self.reference_sequence_repository,
            header.reference_sequences_mut(),
        )?;

        write_header_container(&mut self.inner, &header).await
    }

    /// Writes a SAM header.
    ///
    /// This writes the CRAM magic number, the file definition, and file header using the given SAM
    /// header.
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header).await?;
    ///
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.write_file_definition().await?;
        self.write_file_header(header).await
    }

    /// Writes a CRAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// writer.write_file_definition().await?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header).await?;
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(&header, record).await?;
    ///
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(
        &mut self,
        header: &sam::Header,
        mut record: Record,
    ) -> io::Result<()> {
        use crate::data_container::builder::AddRecordError;

        loop {
            match self.data_container_builder.add_record(record) {
                Ok(_) => {
                    self.record_counter += 1;
                    return Ok(());
                }
                Err(e) => match e {
                    AddRecordError::ContainerFull(r) => {
                        record = r;
                        self.flush(header).await?;
                    }
                    AddRecordError::SliceFull(r) => {
                        record = r;
                    }
                    _ => return Err(io::Error::from(io::ErrorKind::InvalidInput)),
                },
            }
        }
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// writer.write_file_definition().await?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header).await?;
    ///
    /// let record = sam::Record::default();
    /// writer.write_alignment_record(&header, &record).await?;
    ///
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_alignment_record<R>(
        &mut self,
        header: &sam::Header,
        record: &R,
    ) -> io::Result<()>
    where
        R: sam::alignment::Record + ?Sized,
    {
        let record = Record::try_from_alignment_record(header, record)?;
        self.write_record(header, record).await
    }

    async fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::data_container::write_data_container;

        if self.data_container_builder.is_empty() {
            return Ok(());
        }

        let data_container_builder = mem::replace(
            &mut self.data_container_builder,
            DataContainer::builder(self.record_counter),
        );

        let base_count = data_container_builder.base_count();

        let data_container = data_container_builder.build(
            &self.options,
            &self.reference_sequence_repository,
            header,
        )?;

        write_data_container(&mut self.inner, &data_container, base_count).await
    }
}

async fn write_file_definition<W>(
    writer: &mut W,
    file_definition: &FileDefinition,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(&MAGIC_NUMBER).await?;
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
