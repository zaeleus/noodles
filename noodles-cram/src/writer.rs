use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;
use noodles_sam as sam;

use super::{
    container::{self, block, Block},
    num::{write_itf8, write_ltf8},
    MAGIC_NUMBER,
};

// [major, minor]
const FILE_DEFINITION_FORMAT: [u8; 2] = [3, 0];

/// A CRAM writer.
#[derive(Debug)]
pub struct Writer<W> {
    inner: W,
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
        Self { inner }
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

    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        let data = header.to_string().into_bytes();
        let data_len = data.len() as i32;

        let block = Block::new(
            block::CompressionMethod::None,
            block::ContentType::FileHeader as u8,
            0,
            data_len,
            data,
            0,
        );
        let block_len = block.len()? as i32;

        let container_header = container::Header::new(block_len, -1, 0, 0, 0, 0, 0, 1, vec![0], 0);

        write_container_header(&mut self.inner, &container_header)?;
        write_block(&mut self.inner, &block)?;

        Ok(())
    }
}

fn write_container_header<W>(writer: &mut W, header: &container::Header) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length = header.len();
    crc_writer.write_i32::<LittleEndian>(length)?;

    let reference_sequence_id = header.reference_sequence_id();
    write_itf8(&mut crc_writer, reference_sequence_id)?;

    let starting_position_on_the_reference = header.start_position();
    write_itf8(&mut crc_writer, starting_position_on_the_reference)?;

    let alignment_span = header.alignment_span();
    write_itf8(&mut crc_writer, alignment_span)?;

    let number_of_records = header.record_count();
    write_itf8(&mut crc_writer, number_of_records)?;

    let record_counter = header.record_counter();
    write_ltf8(&mut crc_writer, record_counter)?;

    let bases = header.bases();
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

    let content_type = block.content_type();
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
