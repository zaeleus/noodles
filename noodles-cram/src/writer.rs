mod builder;
mod container;
pub(crate) mod data_container;
pub(crate) mod num;
mod options;
pub(crate) mod record;

pub use self::builder::Builder;
pub(crate) use self::options::Options;

use std::{
    io::{self, Write},
    mem,
};

use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::AlignmentRecord;

use self::container::write_container;
use super::{
    container::Container, file_definition::Version, DataContainer, FileDefinition, Record,
    MAGIC_NUMBER,
};

/// A CRAM writer.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_cram as cram;
/// use noodles_fasta as fasta;
/// use noodles_sam as sam;
///
/// let repository = fasta::Repository::default();
/// let header = sam::Header::default();
/// let mut writer = cram::Writer::builder(Vec::new(), repository, &header).build();
/// writer.write_file_definition()?;
///
/// let header = sam::Header::builder().add_comment("noodles-cram").build();
/// writer.write_file_header(&header)?;
///
/// let record = cram::Record::default();
/// writer.write_record(record)?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<'a, W>
where
    W: Write,
{
    inner: W,
    reference_sequence_repository: fasta::Repository,
    header: &'a sam::Header,
    options: Options,
    data_container_builder: crate::data_container::Builder,
    record_counter: i64,
}

impl<'a, W> Writer<'a, W>
where
    W: Write,
{
    /// Creates a CRAM writer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let builder = cram::Writer::builder(Vec::new(), repository, &header);
    /// let writer = builder.build();
    /// ```
    pub fn builder(
        inner: W,
        reference_sequence_repository: fasta::Repository,
        header: &'a sam::Header,
    ) -> Builder<'a, W> {
        Builder::new(inner, reference_sequence_repository, header)
    }

    /// Creates a new CRAM writer with default options.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let writer = cram::Writer::new(Vec::new(), repository, &header);
    /// ```
    pub fn new(
        inner: W,
        reference_sequence_repository: fasta::Repository,
        header: &'a sam::Header,
    ) -> Self {
        Builder::new(inner, reference_sequence_repository, header).build()
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let writer = cram::Writer::new(Vec::new(), repository, &header);
    ///
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Attempts to finish the output stream by writing any pending containers and a final EOF
    /// container.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let mut writer = cram::Writer::new(Vec::new(), repository, &header);
    ///
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.flush()?;
        let eof_container = Container::eof();
        write_container(&mut self.inner, &eof_container)
    }

    /// Writes a CRAM file definition.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let mut writer = cram::Writer::new(Vec::new(), repository, &header);
    ///
    /// writer.write_file_definition()?;
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
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_definition(&mut self) -> io::Result<()> {
        let file_definition = FileDefinition::default();

        // magic number
        self.inner.write_all(MAGIC_NUMBER)?;

        write_format(&mut self.inner, file_definition.version())?;

        self.inner.write_all(file_definition.file_id())?;

        Ok(())
    }

    /// Writes a CRAM file header container.
    ///
    /// The position of the stream is expected to be directly after the file definition.
    ///
    /// Reference sequence dictionary entries must have MD5 checksums (`M5`) set.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let mut writer = cram::Writer::new(Vec::new(), repository, &header);
    ///
    /// writer.write_file_definition()?;
    /// writer.write_file_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        Container::try_from(header)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            .and_then(|container| write_container(&mut self.inner, &container))
    }

    /// Writes a CRAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let repository = fasta::Repository::default();
    /// let header = sam::Header::default();
    /// let mut writer = cram::Writer::new(Vec::new(), repository, &header);
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, mut record: Record) -> io::Result<()> {
        use super::data_container::builder::AddRecordError;

        loop {
            match add_record(
                &mut self.data_container_builder,
                &self.reference_sequence_repository,
                self.header,
                record,
            ) {
                Ok(_) => {
                    self.record_counter += 1;
                    return Ok(());
                }
                Err(e) => match e {
                    AddRecordError::ContainerFull(r) => {
                        record = r;
                        self.flush()?;
                    }
                    AddRecordError::SliceFull(r) => {
                        record = r;
                    }
                },
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
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
            self.header,
        )?;

        let container = Container::try_from_data_container(&data_container, base_count)?;
        write_container(&mut self.inner, &container)?;

        Ok(())
    }
}

impl<'a, W> Drop for Writer<'a, W>
where
    W: Write,
{
    fn drop(&mut self) {
        let _ = self.try_finish();
    }
}

fn add_record(
    data_container_builder: &mut crate::data_container::Builder,
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    record: Record,
) -> Result<(), crate::data_container::builder::AddRecordError> {
    let reference_sequence = record
        .reference_sequence(header.reference_sequences())
        .transpose()
        .expect("invalid reference sequence ID")
        .map(|rs| rs.name())
        .and_then(|name| {
            reference_sequence_repository
                .get(name)
                .transpose()
                .expect("invalid reference sequence")
        })
        .unwrap_or_default();

    data_container_builder.add_record(reference_sequence.as_ref(), record)
}

fn write_format<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let format = [version.major(), version.minor()];
    writer.write_all(&format)
}
