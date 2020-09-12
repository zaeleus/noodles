//! FASTA reader and iterators.

mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Seek, SeekFrom};

use memchr::memchr;

const DEFINITION_PREFIX: u8 = b'>';
const NEWLINE: u8 = b'\n';

/// A FASTA reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a FASTA reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a raw definition line.
    ///
    /// The given buffer will not include the trailing newline. It can subsequently be parsed as a
    /// [`fasta::record::Definition`].
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// definition.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF.
    ///
    /// [`fasta::record::Definition`]: record/definition/struct.Definition.html
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_definition(&mut buf)?;
    ///
    /// assert_eq!(buf, ">sq0");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
    }

    /// Reads a sequence.
    ///
    /// The given buffer consumes a sequence without newlines until another definition or EOF is
    /// reached.
    ///
    /// The position of the stream is expected to be at the start of a sequence, which is directly
    /// after a definition.
    ///
    /// If successful, this returns the number of bytes read from the stream (including the excluded newlines). 
    /// If the number of bytes read is 0, the stream reached EOF (though this case is likely an error).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::Reader::new(&data[..]);
    /// reader.read_definition(&mut String::new())?;
    ///
    /// let mut buf = Vec::new();
    /// reader.read_sequence(&mut buf)?;
    ///
    /// assert_eq!(buf, b"ACGT");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut bytes_read = 0;

        loop {
            let reader_buf = self.inner.fill_buf()?;

            if reader_buf.is_empty() || reader_buf[0] == DEFINITION_PREFIX {
                break;
            }

            let len = match memchr(NEWLINE, reader_buf) {
                Some(i) => {
                    // Discard newline
                    buf.extend(&reader_buf[..i]);
                    i + 1
                }
                None => {
                    buf.extend(reader_buf);
                    reader_buf.len()
                }
            };

            self.inner.consume(len);

            bytes_read += len;
        }

        Ok(bytes_read)
    }

    /// Reads a raw sequence.
    ///
    /// Unlike [read_sequence], the given buffer consumes a sequence including the newlines
    /// until another definition or EOF is reached.
    ///
    /// The position of the stream is expected to be at the start of a sequence, which is directly
    /// after a definition.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF (though this case is likely an error).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\nTGCA\n\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::Reader::new(&data[..]);
    /// reader.read_definition(&mut String::new())?;
    ///
    /// let mut buf = Vec::new();
    /// reader.read_sequence_raw(&mut buf)?;
    ///
    /// assert_eq!(buf, b"ACGT\nTGCA\n\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_sequence_raw(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut bytes_read = 0;

        loop {
            let reader_buf = self.inner.fill_buf()?;

            if reader_buf.is_empty() || reader_buf[0] == DEFINITION_PREFIX {
                break;
            }

            let len = match memchr(NEWLINE, reader_buf) {
                Some(i) => {
                    // Discard newline
                    buf.extend(&reader_buf[..i + 1]);
                    i + 1
                }
                None => {
                    buf.extend(reader_buf);
                    reader_buf.len()
                }
            };

            self.inner.consume(len);

            bytes_read += len;
        }

        Ok(bytes_read)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// definition.
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::Reader::new(&data[..]);
    ///
    /// let mut records = reader.records();
    ///
    /// assert_eq!(records.next().transpose()?, Some(fasta::Record::new(
    ///     fasta::record::Definition::new(String::from("sq0"), None),
    ///     b"ACGT".to_vec(),
    /// )));
    ///
    /// assert_eq!(records.next().transpose()?, Some(fasta::Record::new(
    ///     fasta::record::Definition::new(String::from("sq1"), None),
    ///     b"NNNNNNNNNN".to_vec(),
    /// )));
    ///
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Seek for Reader<R>
where
    R: BufRead + Seek,
{
    /// Seeks the underlying stream to the given position.
    ///
    /// These positions typically come from an associated index, which start at the sequence and
    /// _not_ the definition.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor, Seek, SeekFrom};
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let cursor = Cursor::new(&data[..]);
    /// let mut reader = fasta::Reader::new(cursor);
    /// reader.seek(SeekFrom::Start(14));
    ///
    /// let mut buf = Vec::new();
    /// reader.read_sequence(&mut buf)?;
    ///
    /// assert_eq!(buf, b"NNNNNNNNNN");
    /// # Ok::<(), io::Error>(())
    /// ```
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    let result = reader.read_line(buf);
    buf.pop();
    result
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_read_definition() -> io::Result<()> {
        let data = b">sq0\nACGT\n";
        let mut reader = Reader::new(&data[..]);

        let mut description_buf = String::new();
        reader.read_definition(&mut description_buf)?;

        assert_eq!(description_buf, ">sq0");

        Ok(())
    }

    #[test]
    fn test_read_sequence() -> io::Result<()> {
        let mut sequence_buf = Vec::new();

        let data = b"ACGT\n";
        let mut reader = Reader::new(&data[..]);
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"ACGT\n>sq1\n";
        let mut reader = Reader::new(&data[..]);
        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"NNNN\nNNNN\nNN\n";
        let mut reader = Reader::new(&data[..]);
        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"NNNNNNNNNN");

        Ok(())
    }

    #[test]
    fn test_read_sequence_after_seek() {
        let data = b">sq0\nACGT\n>sq1\nNNNN\n";
        let cursor = Cursor::new(&data[..]);
        let mut reader = Reader::new(cursor);

        reader.seek(SeekFrom::Start(14)).unwrap();

        let mut buf = Vec::new();
        reader.read_sequence(&mut buf).unwrap();

        assert_eq!(buf, b"NNNN");
    }
}
