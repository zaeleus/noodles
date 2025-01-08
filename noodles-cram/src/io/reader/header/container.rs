mod block;
mod header;

use std::io::{self, Read, Take};

use byteorder::{LittleEndian, ReadBytesExt};

use self::block::read_block;
pub(super) use self::header::read_header;

pub(super) struct Reader<R> {
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

    pub(super) fn raw_sam_header_reader(&mut self) -> io::Result<impl Read + '_> {
        let mut reader = read_block(&mut self.inner)?;

        let len = reader.read_i32::<LittleEndian>().and_then(|n| {
            u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        Ok(reader.take(len))
    }
}
