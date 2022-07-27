use std::io;
use std::io::Read;
use super::Index;

/// A GZI index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
    where
        R: Read
{
    /// Creates a GZI index reader.
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a GZI index.
    ///
    /// The position of the [`Read`](std::io::Read) stream is expected to be at the start.
    pub fn read_index(&mut self) -> io::Result<Index> {
        let number_entries = read_little_endian(&mut self.inner)?;
        let mut offsets = Vec::with_capacity(number_entries as usize);

        for _ in 0..number_entries {
            let compressed = read_little_endian(&mut self.inner)?;
            let uncompressed = read_little_endian(&mut self.inner)?;
            offsets.push((compressed, uncompressed));
        }

        if self.inner.read(&mut [0; 1])? > 0 {
            Err(io::Error::new(io::ErrorKind::InvalidData, format!("Trailing data left in read stream, expected {} entries.", number_entries)))
        } else {
            Ok(Index { number_entries, offsets })
        }
    }
}

/// Reads a little endian u64 from the reader.
fn read_little_endian<R>(reader: &mut R) -> io::Result<u64>
    where
      R: Read
{
    let mut buf = [0; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}