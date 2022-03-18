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

use self::container::write_container;
use super::{
    container::Container, file_definition::Version, DataContainer, FileDefinition, Record,
    MAGIC_NUMBER,
};

/// A CRAM writer.
///
/// A call to [`try_finish`] must be made before the writer is dropped.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_cram as cram;
/// use noodles_sam as sam;
///
/// let mut writer = cram::Writer::builder(Vec::new()).build();
/// writer.write_file_definition()?;
///
/// let header = sam::Header::default();
/// writer.write_file_header(&header)?;
///
/// let record = cram::Record::default();
/// writer.write_record(&header, record)?;
///
/// writer.try_finish(&header)?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    reference_sequence_repository: fasta::Repository,
    options: Options,
    data_container_builder: crate::data_container::Builder,
    record_counter: i64,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CRAM writer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let builder = cram::Writer::builder(Vec::new());
    /// let writer = builder.build();
    /// ```
    pub fn builder(inner: W) -> Builder<W> {
        Builder::new(inner)
    }

    /// Creates a new CRAM writer with default options.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::new(inner).build()
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new());
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
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::default();
    /// let mut writer = cram::Writer::new(Vec::new());
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.flush(header)?;
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
    ///
    /// let mut writer = cram::Writer::new(Vec::new());
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
        write_file_definition(&mut self.inner, &file_definition)
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
    /// use noodles_sam as sam;
    ///
    /// let mut writer = cram::Writer::new(Vec::new());
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    ///
    /// writer.try_finish(&header)?;
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
    /// use noodles_sam as sam;
    ///
    ///
    /// let mut writer = cram::Writer::new(Vec::new());
    ///
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(&header, record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &sam::Header, mut record: Record) -> io::Result<()> {
        use super::data_container::builder::AddRecordError;

        loop {
            match self.data_container_builder.add_record(record) {
                Ok(_) => {
                    self.record_counter += 1;
                    return Ok(());
                }
                Err(e) => match e {
                    AddRecordError::ContainerFull(r) => {
                        record = r;
                        self.flush(header)?;
                    }
                    AddRecordError::SliceFull(r) => {
                        record = r;
                    }
                },
            }
        }
    }

    fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
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

        let container = Container::try_from_data_container(&data_container, base_count)?;
        write_container(&mut self.inner, &container)?;

        Ok(())
    }
}

impl<W> sam::AlignmentWriter for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.write_file_definition()?;
        self.write_file_header(header)?;
        Ok(())
    }

    fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::AlignmentRecord,
    ) -> io::Result<()> {
        let r = Record::try_from_alignment_record(header, record)?;
        self.write_record(header, r)
    }

    fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.try_finish(header)
    }
}

fn write_file_definition<W>(writer: &mut W, file_definition: &FileDefinition) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)?;
    write_format(writer, file_definition.version())?;
    writer.write_all(file_definition.file_id())?;
    Ok(())
}

fn write_format<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let format = [version.major(), version.minor()];
    writer.write_all(&format)
}
