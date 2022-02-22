//! BGZF writer.

mod builder;
mod compression_level;

pub use self::{builder::Builder, compression_level::CompressionLevel};

use std::{
    cmp,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::Crc;

use super::{block, gz, BGZF_HEADER_SIZE};

const BGZF_FLG: u8 = 0x04; // FEXTRA
const BGZF_XFL: u8 = 0x00; // none
const BGZF_XLEN: u16 = 6;

const BGZF_SI1: u8 = 0x42;
const BGZF_SI2: u8 = 0x43;
const BGZF_SLEN: u16 = 2;

// § 4.1.2 End-of-file marker (2020-12-03)
pub(crate) static BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, // ID1, ID2
    0x08, // CM = DEFLATE
    0x04, // FLG = FEXTRA
    0x00, 0x00, 0x00, 0x00, // MTIME = 0
    0x00, // XFL = 0
    0xff, // OS = 255 (unknown)
    0x06, 0x00, // XLEN = 6
    0x42, 0x43, // SI1, SI2
    0x02, 0x00, // SLEN = 2
    0x1b, 0x00, // BSIZE = 27
    0x03, 0x00, // CDATA
    0x00, 0x00, 0x00, 0x00, // CRC32 = 0x00000000
    0x00, 0x00, 0x00, 0x00, // ISIZE = 0
];

#[cfg(feature = "libdeflate")]
type CompressionLevelImpl = libdeflater::CompressionLvl;
#[cfg(not(feature = "libdeflate"))]
type CompressionLevelImpl = flate2::Compression;

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
    inner: Option<W>,
    position: u64,
    buf: Vec<u8>,
    compression_level: CompressionLevelImpl,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a BGZF writer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let builder = bgzf::Writer::builder(Vec::new());
    /// let writer = builder.build();
    /// ```
    pub fn builder(inner: W) -> Builder<W> {
        Builder::new(inner)
    }

    /// Creates a writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self::builder(inner).build()
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
        self.inner.as_ref().unwrap()
    }

    fn flush_block(&mut self) -> io::Result<()> {
        let (cdata, crc32, r#isize) = deflate_data(&self.buf, self.compression_level)?;

        let inner = self.inner.as_mut().unwrap();

        write_header(inner, cdata.len())?;
        inner.write_all(&cdata[..])?;
        write_trailer(inner, crc32, r#isize)?;

        let block_size = BGZF_HEADER_SIZE + cdata.len() + gz::TRAILER_SIZE;
        self.position += block_size as u64;

        self.buf.clear();

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

        let inner = self.inner.as_mut().unwrap();
        let result = inner.write_all(BGZF_EOF);

        self.position += BGZF_EOF.len() as u64;

        result
    }

    /// Returns the underlying writer after finishing the output stream.
    ///
    /// This method can only be called once. Any further usage of the writer may result in a panic.
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
        let inner = self.inner.take().unwrap();
        Ok(inner)
    }
}

impl<W> Drop for Writer<W>
where
    W: Write,
{
    fn drop(&mut self) {
        if self.inner.is_some() {
            let _ = self.try_finish();
        }
    }
}

impl<W> Write for Writer<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let max_write_len = cmp::min(
            block::MAX_UNCOMPRESSED_DATA_LENGTH - self.buf.len(),
            buf.len(),
        );

        self.buf.extend_from_slice(&buf[..max_write_len]);

        if self.buf.len() >= block::MAX_UNCOMPRESSED_DATA_LENGTH {
            self.flush()?;
        }

        Ok(max_write_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            Ok(())
        } else {
            self.flush_block()
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

#[cfg(feature = "libdeflate")]
pub(crate) fn deflate_data(
    data: &[u8],
    compression_level: libdeflater::CompressionLvl,
) -> io::Result<(Vec<u8>, u32, u32)> {
    use libdeflater::Compressor;

    let mut encoder = Compressor::new(compression_level);

    let mut compressed_data = Vec::new();

    let max_len = encoder.deflate_compress_bound(data.len());
    compressed_data.resize(max_len, Default::default());

    let len = encoder
        .deflate_compress(data, &mut compressed_data)
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;

    compressed_data.resize(len, Default::default());

    let mut crc = Crc::new();
    crc.update(data);

    Ok((compressed_data, crc.sum(), crc.amount()))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn deflate_data(
    data: &[u8],
    compression_level: flate2::Compression,
) -> io::Result<(Vec<u8>, u32, u32)> {
    use flate2::write::DeflateEncoder;

    let mut encoder = DeflateEncoder::new(Vec::new(), compression_level);
    encoder.write_all(data)?;
    let compressed_data = encoder.finish()?;

    let mut crc = Crc::new();
    crc.update(data);

    Ok((compressed_data, crc.sum(), crc.amount()))
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
