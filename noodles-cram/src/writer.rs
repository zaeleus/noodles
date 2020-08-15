pub mod compression_header;
mod encoding;
pub mod record;
pub mod slice;

use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;
use noodles_sam as sam;

use super::{
    container::{self, block, Block, Container},
    data_container,
    num::{write_itf8, write_ltf8, Itf8},
    DataContainer, Record, MAGIC_NUMBER,
};

// [major, minor]
const FILE_DEFINITION_FORMAT: [u8; 2] = [3, 0];

/// A CRAM writer.
#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    data_container_builder: data_container::Builder,
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
    /// let writer = cram::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            data_container_builder: DataContainer::builder(),
        }
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

    /// Attempts to finish the output stream by writing a final EOF container.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut writer = cram::Writer::new(Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
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
        validate_reference_sequences(header.reference_sequences())?;

        let header_data = header.to_string().into_bytes();
        let header_data_len = header_data.len() as i32;

        let mut data = Vec::new();
        data.write_i32::<LittleEndian>(header_data_len)?;
        data.extend(header_data);

        let block = Block::new(
            block::CompressionMethod::None,
            block::ContentType::FileHeader,
            0,
            data.len() as i32,
            data,
            0,
        );

        let blocks = vec![block];
        let landmarks = vec![0];

        let container_header = container::Header::new(
            0,
            container::ReferenceSequenceId::None,
            0,
            0,
            0,
            0,
            0,
            blocks.len() as Itf8,
            landmarks,
            0,
        );

        let container = Container::new(container_header, blocks);
        self.write_container(&container)?;

        Ok(())
    }

    pub fn write_container(&mut self, container: &Container) -> io::Result<()> {
        let len = container.blocks().iter().map(|b| b.len()).sum();
        write_container_header(&mut self.inner, container.header(), len)?;

        for block in container.blocks() {
            write_block(&mut self.inner, block)?;
        }

        Ok(())
    }

    pub fn write_record(
        &mut self,
        reference_sequence: &[u8],
        mut record: Record,
    ) -> io::Result<()> {
        loop {
            match self
                .data_container_builder
                .add_record(reference_sequence, record)
            {
                Ok(_) => return Ok(()),
                Err(e) => match e {
                    data_container::builder::AddRecordError::ContainerFull(r) => {
                        record = r;

                        let data_container_builder = mem::replace(
                            &mut self.data_container_builder,
                            DataContainer::builder(),
                        );

                        data_container_builder
                            .build()
                            .and_then(|data_container| {
                                Container::try_from_data_container(&data_container)
                            })
                            .and_then(|container| self.write_container(&container))?;
                    }
                },
            }
        }
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
    use noodles_sam::header::reference_sequence::Tag;

    for reference_sequence in reference_sequences.values() {
        if reference_sequence.get(&Tag::Md5Checksum).is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "missing MD5 checksum in reference sequence",
            ));
        }
    }

    Ok(())
}

fn write_container_header<W>(
    writer: &mut W,
    header: &container::Header,
    len: usize,
) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length = len as i32;
    crc_writer.write_i32::<LittleEndian>(length)?;

    let reference_sequence_id = i32::from(header.reference_sequence_id());
    write_itf8(&mut crc_writer, reference_sequence_id)?;

    let starting_position_on_the_reference = header.start_position();
    write_itf8(&mut crc_writer, starting_position_on_the_reference)?;

    let alignment_span = header.alignment_span();
    write_itf8(&mut crc_writer, alignment_span)?;

    let number_of_records = header.record_count();
    write_itf8(&mut crc_writer, number_of_records)?;

    let record_counter = header.record_counter();
    write_ltf8(&mut crc_writer, record_counter)?;

    let bases = header.base_count();
    write_ltf8(&mut crc_writer, bases)?;

    let number_of_blocks = header.block_count();
    write_itf8(&mut crc_writer, number_of_blocks)?;

    let landmarks_len = header.landmarks().len() as i32;
    write_itf8(&mut crc_writer, landmarks_len)?;

    for &pos in header.landmarks() {
        write_itf8(&mut crc_writer, pos)?;
    }

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}

fn write_block<W>(writer: &mut W, block: &Block) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let method = block.compression_method() as u8;
    crc_writer.write_u8(method)?;

    let content_type = block.content_type() as u8;
    crc_writer.write_u8(content_type)?;

    let block_content_id = block.content_id();
    write_itf8(&mut crc_writer, block_content_id)?;

    let size_in_bytes = block.data().len() as i32;
    write_itf8(&mut crc_writer, size_in_bytes)?;

    let uncompressed_data_len = block.uncompressed_len();
    write_itf8(&mut crc_writer, uncompressed_data_len)?;

    crc_writer.write_all(block.data())?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}
