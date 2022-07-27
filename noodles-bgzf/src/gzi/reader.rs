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
            Err(io::Error::new(io::ErrorKind::InvalidData, format!("Trailing data remaining in read stream, expected {} entries.", number_entries)))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_index() -> Result<(), Box<dyn std::error::Error>> {
        let data: &[u8] = &[
            0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // number_entries = 2
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x02, 0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 23810
            0x01, 0x52, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 86529
        ];

        let mut reader = Reader::new(data);

        let actual = reader.read_index()?;
        let expected = Index {
            number_entries: 2, offsets: vec![(4668, 21294), (23810, 86529)]
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_too_many_entries() -> Result<(), Box<dyn std::error::Error>> {
        let data: &[u8] = &[
            0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // number_entries = 3
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x02, 0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 23810
            0x01, 0x52, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 86529
        ];

        let mut reader = Reader::new(data);

        let actual = reader.read_index().map_err(|err| err.kind());
        let expected = Err(io::ErrorKind::UnexpectedEof);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_trailing_data() -> Result<(), Box<dyn std::error::Error>> {
        let data: &[u8] = &[
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // number_entries = 1
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x00
        ];

        let mut reader = Reader::new(data);

        let actual = reader.read_index().map_err(|err| err.kind());
        let expected = Err(io::ErrorKind::InvalidData);

        assert_eq!(actual, expected);

        Ok(())
    }
}
