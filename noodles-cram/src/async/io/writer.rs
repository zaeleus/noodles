//! Async CRAM writer.

mod builder;
mod container;
mod header;

use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;
use self::{
    container::write_container,
    header::{write_file_definition, write_file_header, write_header},
};
use crate::{
    FileDefinition,
    io::writer::{Options, Record},
};

/// An async CRAM writer.
///
/// A call to [`Self::shutdown`] must be made before the writer is dropped.
pub struct Writer<W> {
    inner: W,
    reference_sequence_repository: fasta::Repository,
    options: Options,
    records: Vec<Record>,
    record_counter: u64,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
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

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
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

    /// Attempts to shutdown the output stream by writing any pending containers and a final EOF
    /// container.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
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
        write_eof_container(&mut self.inner, self.options.version).await
    }

    /// Writes a CRAM file definition.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let mut writer = cram::r#async::io::Writer::new(io::sink());
    /// writer.write_file_definition().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_file_definition(&mut self) -> io::Result<()> {
        let file_definition = FileDefinition::new(self.options.version, Default::default());
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
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        write_file_header(
            &mut self.inner,
            &self.reference_sequence_repository,
            header,
            self.options.version,
            self.options.reference_required,
        )
        .await
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
        let file_definition = FileDefinition::new(self.options.version, Default::default());

        write_header(
            &mut self.inner,
            &self.reference_sequence_repository,
            &file_definition,
            header,
            self.options.reference_required,
        )
        .await
    }

    /// Writes a CRAM record.
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
    /// let record = cram::Record::default();
    /// writer.write_record(&header, &record).await?;
    ///
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(
        &mut self,
        header: &sam::Header,
        record: &crate::Record<'_>,
    ) -> io::Result<()> {
        let writer_record = Record::from_cram_record(record, self.options.strip_md_nm)?;
        self.add_record(header, writer_record).await
    }

    async fn add_record(&mut self, header: &sam::Header, record: Record) -> io::Result<()> {
        self.records.push(record);

        if self.records.len() >= self.records.capacity() {
            self.flush(header).await?;
        }

        Ok(())
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
    pub async fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        let record = Record::try_from_alignment_record_with_options(
            header,
            record,
            self.options.strip_md_nm,
        )?;
        self.add_record(header, record).await
    }

    async fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        write_container(
            &mut self.inner,
            &self.reference_sequence_repository,
            &self.options,
            header,
            self.record_counter,
            &mut self.records,
        )
        .await?;

        let record_count = u64::try_from(self.records.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.record_counter += record_count;

        self.records.clear();

        Ok(())
    }
}
