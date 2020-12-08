use std::{
    cmp,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::{write::DeflateEncoder, Compression, Crc};

use super::{block, gz, BGZF_HEADER_SIZE};

const MAX_COMPRESSED_DATA_LENGTH: usize = block::MAX_LENGTH + 1; // bytes

const BGZF_FLG: u8 = 0x04; // FEXTRA
const BGZF_XFL: u8 = 0x00; // none
const BGZF_XLEN: u16 = 6;

const BGZF_SI1: u8 = 0x42;
const BGZF_SI2: u8 = 0x43;
const BGZF_SLEN: u16 = 2;

// Sequence Alignment/Map Format Specification § 4.1.2 (accessed 2020-04-15)
static BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// A BZGF writer.
///
/// This implements [`std::io::Write`], consuming uncompressed data and emitting compressed data.
///
/// # Examples
///
/// ```
/// # use std::io::{self, Write};
/// use noodles_bgzf as bgzf;
///
/// let mut writer = bgzf::Writer::new(Vec::new());
/// writer.write_all(b"noodles-bgzf")?;
///
/// let data = writer.finish()?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    encoder: DeflateEncoder<Vec<u8>>,
    crc: Crc,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            encoder: DeflateEncoder::new(Vec::new(), Compression::default()),
            crc: Crc::new(),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    fn flush_block(&mut self) -> io::Result<()> {
        self.encoder.try_finish()?;
        let data = self.encoder.get_ref();

        write_header(&mut self.inner, data.len())?;
        self.inner.write_all(&data[..])?;
        write_trailer(&mut self.inner, self.crc.sum(), self.crc.amount())?;

        self.encoder.reset(Vec::new())?;
        self.crc.reset();

        Ok(())
    }

    /// Attempts to finish the output stream by flushing any remaining buffers.
    ///
    /// This then appends the final BGZF EOF block.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut writer = bgzf::Writer::new(Vec::new());
    /// writer.write_all(b"noodles-bgzf")?;
    ///
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.flush()?;
        self.inner.write_all(BGZF_EOF)
    }

    /// Returns the underlying writer after finishing the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut writer = bgzf::Writer::new(Vec::new());
    /// writer.write_all(b"noodles-bgzf")?;
    ///
    /// let data = writer.finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn finish(mut self) -> io::Result<W> {
        self.try_finish()?;
        Ok(self.inner)
    }
}

impl<W> Write for Writer<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let total_uncompressed_bytes_written = self.crc.amount() as usize;

        // Only the uncompressed size is tracked, and the assumption is that the uncompressed size
        // will always be less than the compressed size.
        if total_uncompressed_bytes_written >= MAX_COMPRESSED_DATA_LENGTH {
            self.flush()?;
            return Err(io::Error::from(io::ErrorKind::Interrupted));
        }

        let bytes_to_be_written = cmp::min(
            (MAX_COMPRESSED_DATA_LENGTH - total_uncompressed_bytes_written) as usize,
            buf.len(),
        );
        let bytes_written = self.encoder.write(&buf[..bytes_to_be_written])?;
        self.crc.update(&buf[..bytes_written]);

        Ok(bytes_written)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.crc.amount() > 0 {
            self.flush_block()
        } else {
            Ok(())
        }
    }
}

fn write_header<W>(writer: &mut W, cdata_len: usize) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&gz::MAGIC_NUMBER)?;
    writer.write_u8(gz::CompressionMethod::Deflate as u8)?;
    writer.write_u8(BGZF_FLG)?;
    writer.write_u32::<LittleEndian>(gz::MTIME_NONE)?;
    writer.write_u8(BGZF_XFL)?;
    writer.write_u8(gz::OperatingSystem::Unknown as u8)?;
    writer.write_u16::<LittleEndian>(BGZF_XLEN)?;

    writer.write_u8(BGZF_SI1)?;
    writer.write_u8(BGZF_SI2)?;
    writer.write_u16::<LittleEndian>(BGZF_SLEN)?;

    let bsize = (cdata_len + BGZF_HEADER_SIZE + gz::TRAILER_SIZE - 1) as u16;
    writer.write_u16::<LittleEndian>(bsize)?;

    Ok(())
}

fn write_trailer<W>(writer: &mut W, checksum: u32, uncompressed_size: u32) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(checksum)?;
    writer.write_u32::<LittleEndian>(uncompressed_size)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_finish() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_all(b"noodles")?;

        let data = writer.finish()?;
        let eof_start = data.len() - BGZF_EOF.len();

        assert_eq!(&data[eof_start..], BGZF_EOF);

        Ok(())
    }
}
