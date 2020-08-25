mod block;
pub mod compression_header;
mod container;
mod encoding;
pub mod record;
pub mod slice;

use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{
    container::{Block, Container},
    data_container,
    num::Itf8,
    DataContainer, Record, MAGIC_NUMBER,
};

use self::block::write_block;

// [major, minor]
const FILE_DEFINITION_FORMAT: [u8; 2] = [3, 0];

const RECORD_COUNTER_START: i64 = 0;

/// A CRAM writer.
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

    /// Writes a CRAM file defintion.
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
        // magic number
        self.inner.write_all(MAGIC_NUMBER)?;

        // format (major, minor)
        self.inner.write_all(&FILE_DEFINITION_FORMAT)?;

        // File ID is currently blank.
        let file_id = [0; 20];
        self.inner.write_all(&file_id)?;

        Ok(())
    }

    /// Writes a CRAM file header container.
    ///
    /// The position of the stream is expected to be directly after the file definition.
    ///
    /// Reference sequence dictionary entries must have MD5 checksums (`M5`) set.
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        use crate::container::block::ContentType;

        validate_reference_sequences(header.reference_sequences())?;

        let header_data = header.to_string().into_bytes();
        let header_data_len = header_data.len() as i32;

        let mut data = Vec::new();
        data.write_i32::<LittleEndian>(header_data_len)?;
        data.extend(header_data);

        let block = Block::builder()
            .set_content_type(ContentType::FileHeader)
            .set_uncompressed_len(data.len() as Itf8)
            .set_data(data)
            .build();

        let blocks = vec![block];
        let landmarks = vec![0];

        let len = blocks.iter().map(|b| b.len() as i32).sum();

        let container_header = crate::container::Header::builder()
            .set_length(len)
            .set_reference_sequence_id(crate::container::ReferenceSequenceId::None)
            .set_block_count(blocks.len() as Itf8)
            .set_landmarks(landmarks)
            .build();

        let container = Container::new(container_header, blocks);
        self.write_container(&container)?;

        Ok(())
    }

    pub fn write_container(&mut self, container: &Container) -> io::Result<()> {
        self::container::write_header(&mut self.inner, container.header())?;

        for block in container.blocks() {
            write_block(&mut self.inner, block)?;
        }

        Ok(())
    }

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

        data_container_builder
            .build(&self.reference_sequences)
            .and_then(|data_container| Container::try_from_data_container(&data_container))
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

fn validate_reference_sequences(
    reference_sequences: &sam::header::ReferenceSequences,
) -> io::Result<()> {
    for reference_sequence in reference_sequences.values() {
        if reference_sequence.md5_checksum().is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "missing MD5 checksum in reference sequence",
            ));
        }
    }

    Ok(())
}

fn add_record(
    data_container_builder: &mut data_container::Builder,
    reference_sequences: &[fasta::Record],
    record: Record,
) -> Result<(), data_container::builder::AddRecordError> {
    let reference_sequence = record
        .reference_sequence_id()
        .and_then(|id| reference_sequences.get(id as usize))
        .map(|rs| rs.sequence())
        .unwrap_or_default();

    data_container_builder.add_record(reference_sequence, record)
}
