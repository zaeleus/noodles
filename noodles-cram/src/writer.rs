mod block;
pub(crate) mod compression_header;
mod container;
mod encoding;
pub(crate) mod num;
pub(crate) mod record;
pub(crate) mod slice;

use std::{
    convert::TryFrom,
    io::{self, Write},
    mem,
};

use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{
    container::Container, data_container, file_definition::Version, DataContainer, FileDefinition,
    Record, MAGIC_NUMBER,
};

use self::block::write_block;

const RECORD_COUNTER_START: i64 = 0;

/// A CRAM writer.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_cram as cram;
/// use noodles_sam as sam;
///
/// let mut writer = cram::Writer::new(Vec::new(), Vec::new());
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
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    reference_sequences: Vec<fasta::Record>,
    data_container_builder: data_container::Builder,
    record_counter: i64,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a new CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new(), Vec::new());
    /// ```
    pub fn new(inner: W, reference_sequences: Vec<fasta::Record>) -> Self {
        Self {
            inner,
            reference_sequences,
            data_container_builder: DataContainer::builder(RECORD_COUNTER_START),
            record_counter: RECORD_COUNTER_START,
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new(), Vec::new());
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
    /// let mut writer = cram::Writer::new(Vec::new(), Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.flush()?;
        let eof_container = Container::eof();
        self.write_container(&eof_container)
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
    /// let mut writer = cram::Writer::new(Vec::new(), Vec::new());
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
    /// use noodles_sam as sam;
    ///
    /// let mut writer = cram::Writer::new(Vec::new(), Vec::new());
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        Container::try_from(header)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            .and_then(|container| self.write_container(&container))
    }

    /// Writes a CRAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut writer = cram::Writer::new(Vec::new(), Vec::new());
    /// let record = cram::Record::default();
    /// writer.write_record(record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, mut record: Record) -> io::Result<()> {
        loop {
            match add_record(
                &mut self.data_container_builder,
                &self.reference_sequences,
                record,
            ) {
                Ok(_) => {
                    self.record_counter += 1;
                    return Ok(());
                }
                Err(e) => match e {
                    data_container::builder::AddRecordError::ContainerFull(r) => {
                        record = r;
                        self.flush()?;
                    }
                    data_container::builder::AddRecordError::SliceFull(r) => {
                        record = r;
                    }
                },
            }
        }
    }

    fn write_container(&mut self, container: &Container) -> io::Result<()> {
        self::container::write_header(&mut self.inner, container.header())?;

        for block in container.blocks() {
            write_block(&mut self.inner, block)?;
        }

        Ok(())
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

        data_container_builder
            .build(&self.reference_sequences)
            .and_then(|data_container| {
                Container::try_from_data_container(&data_container, base_count)
            })
            .and_then(|container| self.write_container(&container))
    }
}

impl<W> Drop for Writer<W>
where
    W: Write,
{
    fn drop(&mut self) {
        let _ = self.try_finish();
    }
}

fn add_record(
    data_container_builder: &mut data_container::Builder,
    reference_sequences: &[fasta::Record],
    record: Record,
) -> Result<(), data_container::builder::AddRecordError> {
    let reference_sequence = record
        .reference_sequence_id()
        .map(i32::from)
        .and_then(|id| reference_sequences.get(id as usize))
        .map(|rs| rs.sequence())
        .unwrap_or_default();

    data_container_builder.add_record(reference_sequence, record)
}

fn write_format<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let format = [version.major(), version.minor()];
    writer.write_all(&format)
}
