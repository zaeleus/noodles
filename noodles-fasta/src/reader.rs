//! FASTA reader and iterators.

mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use memchr::memchr;
use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Region};

use super::{fai, Record};

pub(crate) const DEFINITION_PREFIX: u8 = b'>';
pub(crate) const NEWLINE: u8 = b'\n';

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

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
    /// [`crate::record::Definition`].
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// definition.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF.
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
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF (though this case is likely an error).
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
                    let line = &reader_buf[..i];

                    if line.ends_with(&[CARRIAGE_RETURN as u8]) {
                        let end = line.len() - 1;
                        buf.extend(&line[..end]);
                    } else {
                        buf.extend(line);
                    }

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

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// use noodles_fasta as fasta;
    ///
    /// let data = Vec::new();
    /// let reader = fasta::Reader::new(bgzf::Reader::new(&data[..]));
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(virtual_position.compressed(), 0);
    /// assert_eq!(virtual_position.uncompressed(), 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<R>
where
    R: BufRead + Seek,
{
    /// Returns a record of the given region.
    ///
    /// This currently only supports querying the entirety of a reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_core::Region;
    /// use noodles_fasta::{self as fasta, fai};
    ///
    /// let data = b">sq0\nNNNN\n>sq1\nACGT\n>sq2\nNNNN\n";
    /// let index = vec![
    ///     fai::Record::new(String::from("sq0"), 4, 5, 4, 5),
    ///     fai::Record::new(String::from("sq1"), 4, 15, 4, 5),
    ///     fai::Record::new(String::from("sq2"), 4, 25, 4, 5),
    /// ];
    ///
    /// let mut reader = fasta::Reader::new(Cursor::new(data));
    ///
    /// let region = Region::mapped("sq1", ..);
    /// let record = reader.query(&index, &region)?;
    ///
    /// assert_eq!(record, fasta::Record::new(
    ///     fasta::record::Definition::new(String::from("sq1"), None),
    ///     b"ACGT".to_vec(),
    /// ));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn query(&mut self, index: &[fai::Record], region: &Region) -> io::Result<Record> {
        use crate::record::Definition;

        let (i, _) = resolve_region(index, region)?;
        let index_record = &index[i];

        let pos = index_record.offset();
        self.seek(SeekFrom::Start(pos))?;

        let definition = Definition::new(region.to_string(), None);
        let mut sequence = Vec::new();

        self.read_sequence(&mut sequence)?;

        Ok(Record::new(definition, sequence))
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read + Seek,
{
    /// Seeks the underlying BGZF stream to the given virtual position.
    ///
    /// Virtual positions typically come from an associated index.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bgzf as bgzf;
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.fa.gz")
    ///     .map(bgzf::Reader::new)
    ///     .map(fasta::Reader::new)?;
    ///
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }
}

impl<R> Seek for Reader<R>
where
    R: Read + Seek,
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

// Reads all bytes until a line feed ('\n') or EOF is reached.
//
// The buffer will not include the trailing newline ('\n' or '\r\n').
pub(crate) fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_line(buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

fn resolve_region(index: &[fai::Record], region: &Region) -> io::Result<(usize, Interval)> {
    if let Some(r) = region.as_mapped() {
        let i = index
            .iter()
            .position(|r| r.reference_sequence_name() == region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", r.name()),
                )
            })?;

        Ok((i, r.interval()))
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "region is not mapped",
        ))
    }
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
    fn test_read_sequence_with_crlf() -> io::Result<()> {
        let mut sequence_buf = Vec::new();

        let data = b"ACGT\r\n";
        let mut reader = Reader::new(&data[..]);
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"ACGT\r\n>sq1\r\n";
        let mut reader = Reader::new(&data[..]);
        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"NNNN\r\nNNNN\r\nNN\r\n";
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

    #[test]
    fn test_read_line() -> io::Result<()> {
        let mut buf = String::new();

        let data = b"noodles\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        let data = b"noodles\r\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        let data = b"noodles";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        Ok(())
    }
}
