//! FASTA indexer.

use std::{
    error::Error,
    fmt,
    io::{self, BufRead},
};

use memchr::memchr;

use super::reader::{read_line, DEFINITION_PREFIX};
use crate::{
    fai::Record,
    record::definition::{Definition, ParseError},
};

/// A FASTA indexer.
pub struct Indexer<R> {
    inner: R,
    offset: u64,
}

impl<R> Indexer<R>
where
    R: BufRead,
{
    /// Creates a FASTA indexer.
    pub fn new(inner: R) -> Self {
        Self { inner, offset: 0 }
    }

    /// Consumes a single sequence line.
    ///
    /// If successful, this returns the number of bytes read from the stream (i.e., the line width)
    /// and the number of bases in the line. If the number of bytes read is 0, the entire sequence
    /// of the current record was read.
    fn consume_sequence_line(&mut self) -> io::Result<(usize, usize)> {
        consume_sequence_line(&mut self.inner)
    }

    /// Indexes a raw FASTA record.
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// definition.
    ///
    /// # Errors
    ///
    /// An error is returned if the record fails to be completely read. This includes when
    ///
    ///   * the stream is not at the start of a definition;
    ///   * the record is missing a sequence;
    ///   * the sequence lines have a different number of bases, excluding the last line;
    ///   * or the sequence lines are not the same length, excluding the last line.
    pub fn index_record(&mut self) -> Result<Option<Record>, IndexError> {
        let definition = match self.read_definition() {
            Ok(None) => return Ok(None),
            Ok(Some(d)) => d,
            Err(e) => return Err(e.into()),
        };

        let offset = self.offset;
        let mut length = 0;

        let (line_width, line_bases) = self.consume_sequence_line()?;
        let (mut prev_line_width, mut prev_line_bases) = (line_width, line_bases);

        loop {
            self.offset += prev_line_width as u64;
            length += prev_line_bases;

            match self.consume_sequence_line() {
                Ok((0, _)) => break,
                Ok((bytes_read, base_count)) => {
                    if line_bases != prev_line_bases {
                        return Err(IndexError::InvalidLineBases(line_bases, prev_line_bases));
                    } else if line_width != prev_line_width {
                        return Err(IndexError::InvalidLineWidth(line_width, prev_line_width));
                    }

                    prev_line_width = bytes_read;
                    prev_line_bases = base_count;
                }
                Err(e) => return Err(IndexError::IoError(e)),
            }
        }

        if length == 0 {
            return Err(IndexError::EmptySequence(self.offset));
        }

        let record = Record::new(
            definition.name(),
            length as u64,
            offset,
            line_bases as u64,
            line_width as u64,
        );

        Ok(Some(record))
    }

    fn read_definition(&mut self) -> io::Result<Option<Definition>> {
        let mut buf = String::new();

        match read_line(&mut self.inner, &mut buf) {
            Ok(0) => return Ok(None),
            Ok(n) => self.offset += n as u64,
            Err(e) => return Err(e),
        }

        buf.parse()
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

fn consume_sequence_line<R>(reader: &mut R) -> io::Result<(usize, usize)>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    fn count_bases(buf: &[u8]) -> usize {
        if buf.ends_with(&[CARRIAGE_RETURN]) {
            buf.len() - 1
        } else {
            buf.len()
        }
    }

    let mut bytes_read = 0;
    let mut base_count = 0;
    let mut is_eol = false;

    loop {
        let src = reader.fill_buf()?;

        if is_eol || src.is_empty() || src[0] == DEFINITION_PREFIX {
            break;
        }

        let (chunk_len, chunk_base_count) = match memchr(LINE_FEED, src) {
            Some(i) => {
                is_eol = true;
                (i + 1, count_bases(&src[..i]))
            }
            None => (src.len(), count_bases(src)),
        };

        reader.consume(chunk_len);

        bytes_read += chunk_len;
        base_count += chunk_base_count;
    }

    Ok((bytes_read, base_count))
}

#[derive(Debug)]
pub enum IndexError {
    EmptySequence(u64),
    InvalidDefinition(ParseError),
    InvalidLineBases(usize, usize),
    InvalidLineWidth(usize, usize),
    IoError(io::Error),
}

impl Error for IndexError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::EmptySequence(_) => None,
            Self::InvalidDefinition(e) => Some(e),
            Self::InvalidLineBases(..) => None,
            Self::InvalidLineWidth(..) => None,
            Self::IoError(e) => Some(e),
        }
    }
}

impl fmt::Display for IndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequence(offset) => write!(f, "empty sequence at offset {offset}"),
            Self::InvalidDefinition(e) => e.fmt(f),
            Self::InvalidLineBases(expected, actual) => {
                write!(f, "invalid line bases: expected {expected}, got {actual}")
            }
            Self::InvalidLineWidth(expected, actual) => {
                write!(f, "invalid line width: expected {expected}, got {actual}")
            }
            Self::IoError(e) => e.fmt(f),
        }
    }
}

impl From<io::Error> for IndexError {
    fn from(error: io::Error) -> Self {
        Self::IoError(error)
    }
}

impl From<ParseError> for IndexError {
    fn from(error: ParseError) -> Self {
        Self::InvalidDefinition(error)
    }
}

impl From<IndexError> for io::Error {
    fn from(error: IndexError) -> Self {
        match error {
            IndexError::IoError(e) => e,
            _ => Self::new(io::ErrorKind::InvalidInput, error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_record() -> Result<(), IndexError> {
        let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
        let mut indexer = Indexer::new(&data[..]);

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new("sq0", 4, 5, 4, 5)));

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new("sq1", 10, 15, 4, 5)));

        assert!(indexer.index_record()?.is_none());

        Ok(())
    }

    #[test]
    fn test_index_record_with_invalid_line_bases() {
        let data = b">sq0\nACGT\nACG\nACGT\nAC\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineBases(4, 3))
        ));
    }

    #[test]
    fn test_index_record_with_invalid_line_width() {
        let data = b">sq0\nACGT\nACGT\r\nACGT\nAC\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineWidth(5, 6))
        ));
    }

    #[test]
    fn test_index_record_with_empty_sequence() {
        let data = b">sq0\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::EmptySequence(5))
        ));
    }

    #[test]
    fn test_consume_sequence_line() -> io::Result<()> {
        use std::io::BufReader;

        let data = b"ACGT\nNNNN\n";
        let mut reader = &data[..];
        let (len, base_count) = consume_sequence_line(&mut reader)?;
        assert_eq!(len, 5);
        assert_eq!(base_count, 4);

        let data = b"ACGT\r\nNNNN\r\n";
        let mut reader = &data[..];
        let (len, base_count) = consume_sequence_line(&mut reader)?;
        assert_eq!(len, 6);
        assert_eq!(base_count, 4);

        let data = b"ACGT\r\nNNNN\r\n";
        let mut reader = BufReader::with_capacity(3, &data[..]);
        let (len, base_count) = consume_sequence_line(&mut reader)?;
        assert_eq!(len, 6);
        assert_eq!(base_count, 4);

        Ok(())
    }
}
