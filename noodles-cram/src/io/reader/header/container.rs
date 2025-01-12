//! CRAM header container reader.

mod block;
mod header;
pub mod sam_header;

use std::io::{self, Read, Take};

use byteorder::{LittleEndian, ReadBytesExt};

use self::block::read_block;
pub(super) use self::header::read_header;

/// A CRAM header container reader.
pub struct Reader<R> {
    inner: Take<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R, len: u64) -> Self {
        Self {
            inner: inner.take(len),
        }
    }

    /// Returns a raw SAM header reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header text, e.g., using
    /// [`sam_header::Reader::discard_to_end`].
    pub fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<impl Read + '_>> {
        let mut reader = read_block(&mut self.inner)?;

        let len = reader.read_i32::<LittleEndian>().and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        Ok(sam_header::Reader::new(reader, len))
    }

    /// Discards all input until EOF.
    pub fn discard_to_end(&mut self) -> io::Result<u64> {
        io::copy(&mut self.inner, &mut io::sink())
    }
}
