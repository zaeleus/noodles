//! CRAM header container reader.

mod block;
mod header;
pub mod sam_header;

use std::io::{self, Read, Take};

use self::block::read_block;
pub(super) use self::header::read_header;
use crate::{file_definition::Version, io::reader::num::read_i32_le};

/// A CRAM header container reader.
pub struct Reader<R> {
    inner: Take<R>,
    version: Version,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R, len: u64, version: Version) -> Self {
        Self {
            inner: inner.take(len),
            version,
        }
    }

    /// Returns a raw SAM header reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header text, e.g., using
    /// [`sam_header::Reader::discard_to_end`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{io::Read, fs::File};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    /// header_reader.read_file_id()?;
    ///
    /// let mut container_reader = header_reader.container_reader()?;
    ///
    /// let buf = {
    ///     let mut buf = Vec::new();
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    ///     raw_sam_header_reader.read_to_end(&mut buf)?;
    ///     raw_sam_header_reader.discard_to_end()?;
    ///     buf
    /// };
    ///
    /// container_reader.discard_to_end()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<impl Read + '_>> {
        let mut reader = read_block(&mut self.inner, self.version)?;

        let len = read_i32_le(&mut reader).and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        Ok(sam_header::Reader::new(reader, len))
    }

    /// Discards all input until EOF.
    pub fn discard_to_end(&mut self) -> io::Result<u64> {
        io::copy(&mut self.inner, &mut io::sink())
    }
}
