//! BGZF writer.

mod builder;
mod compression_level;
mod frame;

pub use self::{builder::Builder, compression_level::CompressionLevel};

use std::{
    cmp,
    io::{self, Write},
};

use super::{gz, VirtualPosition, BGZF_HEADER_SIZE, BGZF_MAX_ISIZE};

// The max DEFLATE overhead for 65536 bytes of data at compression level 0.
//
// For zlib (and derivatives) and libdeflate, this is 10 bytes; and for miniz_oxide, 15 bytes.
const COMPRESSION_LEVEL_0_OVERHEAD: usize = 15;

// The max size of the write buffer.
//
// The buffer that uses this size is the uncompressed data that is staged to be written as a BGZF
// block. It is slightly smaller than the max allowed ISIZE to compensate for the gzip format and
// DEFLATE overheads.
pub(crate) const MAX_BUF_SIZE: usize =
    BGZF_MAX_ISIZE - BGZF_HEADER_SIZE - gz::TRAILER_SIZE - COMPRESSION_LEVEL_0_OVERHEAD;

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
pub(crate) type CompressionLevelImpl = libdeflater::CompressionLvl;
#[cfg(not(feature = "libdeflate"))]
pub(crate) type CompressionLevelImpl = flate2::Compression;

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
    /// Creates a writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_with_writer(inner)
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

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(mut self) -> W {
        self.inner.take().unwrap()
    }

    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// assert_eq!(writer.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Panics
    ///
    /// This panics if the stream flushed >= 256 TiB of compressed data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::Writer::new(Vec::new());
    /// assert_eq!(writer.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        // SAFETY: The uncompressed buffer is guaranteed to be <= `MAX_UNCOMPRESSED_POSITION`.
        let uncompressed_position = self.buf.len() as u16;
        VirtualPosition::try_from((self.position, uncompressed_position)).unwrap()
    }

    fn flush_block(&mut self) -> io::Result<()> {
        use self::frame::{write_header, write_trailer};
        use crate::deflate;

        let (cdata, crc32, r#isize) = deflate::encode(&self.buf, self.compression_level)?;

        let inner = self.inner.as_mut().unwrap();
        let block_size = BGZF_HEADER_SIZE + cdata.len() + gz::TRAILER_SIZE;

        write_header(inner, block_size)?;
        inner.write_all(&cdata[..])?;
        write_trailer(inner, crc32, r#isize)?;

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
        let max_write_len = cmp::min(MAX_BUF_SIZE - self.buf.len(), buf.len());

        self.buf.extend_from_slice(&buf[..max_write_len]);

        if self.buf.len() >= MAX_BUF_SIZE {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_virtual_position() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        assert_eq!(writer.virtual_position(), VirtualPosition::from(0));

        writer.write_all(b"noodles")?;

        assert_eq!(
            writer.virtual_position(),
            VirtualPosition::try_from((0, 7))?
        );

        writer.flush()?;

        assert_eq!(
            writer.virtual_position(),
            VirtualPosition::try_from((writer.get_ref().len() as u64, 0))?
        );

        Ok(())
    }

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
