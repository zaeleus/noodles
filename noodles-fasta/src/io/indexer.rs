//! FASTA indexer.

use std::{
    error::Error,
    fmt,
    io::{self, BufRead},
    num::NonZero,
};

use memchr::memchr;

use super::reader::{DEFINITION_PREFIX, read_line};
use crate::{fai::Record, record::definition::Definition};

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
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    /// let mut indexer = fasta::io::Indexer::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner, offset: 0 }
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
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use noodles_fasta::{self as fasta, fai};
    ///
    /// let src = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut indexer = fasta::io::Indexer::new(&src[..]);
    ///
    /// let mut records = Vec::new();
    ///
    /// while let Some(record) = indexer.index_record()? {
    ///     records.push(record);
    /// }
    ///
    /// let line_base_count = const { NonZero::new(4).unwrap() };
    /// let line_width = const { NonZero::new(5).unwrap() };
    /// let expected = [
    ///     fai::Record::new("sq0", 4, 5, line_base_count, line_width),
    ///     fai::Record::new("sq1", 10, 15, line_base_count, line_width),
    /// ];
    ///
    /// assert_eq!(records, expected);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn index_record(&mut self) -> Result<Option<Record>, IndexError> {
        let Some(definition) = self.read_definition()? else {
            return Ok(None);
        };

        let offset = self.offset;

        let (expected_line_width, expected_line_base_count) = self.consume_sequence_line()?;
        let mut base_count = expected_line_base_count;

        if base_count == 0 {
            return Err(IndexError::EmptySequence(self.offset));
        }

        loop {
            let (line_width, line_base_count) = self.consume_sequence_line()?;

            base_count += line_base_count;

            let is_valid_last_sequence_line = is_last_sequence_line(&mut self.inner)?
                && line_width <= expected_line_width
                && line_base_count <= expected_line_base_count;

            if is_valid_last_sequence_line {
                break;
            }

            if line_base_count != expected_line_base_count {
                return Err(IndexError::InvalidLineBases(
                    line_base_count,
                    expected_line_base_count,
                ));
            } else if line_width != expected_line_width {
                return Err(IndexError::InvalidLineWidth(
                    line_width,
                    expected_line_width,
                ));
            }
        }

        let line_base_count = u64::try_from(expected_line_base_count)
            .and_then(NonZero::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let line_width = u64::try_from(expected_line_width)
            .and_then(NonZero::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let record = Record::new(
            definition.name(),
            base_count as u64,
            offset,
            line_base_count,
            line_width,
        );

        Ok(Some(record))
    }

    fn read_definition(&mut self) -> io::Result<Option<Definition>> {
        use super::reader::definition::parse_definition;

        let mut buf = String::new();

        match read_line(&mut self.inner, &mut buf) {
            Ok(0) => return Ok(None),
            Ok(n) => self.offset += n as u64,
            Err(e) => return Err(e),
        }

        let (name, description) = parse_definition(buf.as_bytes())?;

        let description = if description.is_empty() {
            None
        } else {
            Some(description.into())
        };

        Ok(Some(Definition::new(name, description)))
    }

    fn consume_sequence_line(&mut self) -> io::Result<(usize, usize)> {
        let (line_width, line_base_count) = consume_sequence_line(&mut self.inner)?;
        self.offset += line_width as u64;
        Ok((line_width, line_base_count))
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

fn is_last_sequence_line<R>(reader: &mut R) -> io::Result<bool>
where
    R: BufRead,
{
    let src = reader.fill_buf()?;
    Ok(src.is_empty() || src[0] == DEFINITION_PREFIX)
}

#[derive(Debug)]
pub enum IndexError {
    Io(io::Error),
    EmptySequence(u64),
    InvalidLineBases(usize, usize),
    InvalidLineWidth(usize, usize),
}

impl Error for IndexError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::EmptySequence(_) => None,
            Self::InvalidLineBases(..) => None,
            Self::InvalidLineWidth(..) => None,
        }
    }
}

impl fmt::Display for IndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => e.fmt(f),
            Self::EmptySequence(offset) => write!(f, "empty sequence at offset {offset}"),
            Self::InvalidLineBases(actual, expected) => {
                write!(f, "invalid line bases: expected {expected}, got {actual}")
            }
            Self::InvalidLineWidth(actual, expected) => {
                write!(f, "invalid line width: expected {expected}, got {actual}")
            }
        }
    }
}

impl From<io::Error> for IndexError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<IndexError> for io::Error {
    fn from(error: IndexError) -> Self {
        match error {
            IndexError::Io(e) => e,
            _ => Self::new(io::ErrorKind::InvalidInput, error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_record_with_invalid_line_bases() {
        let data = b">sq0\nACGT\nACG\nACGT\nAC\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineBases(3, 4))
        ));

        let data = b">sq0\nACGT\nACGTN\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineBases(5, 4))
        ));
    }

    #[test]
    fn test_index_record_with_invalid_line_width() {
        let data = b">sq0\nACGT\nACGT\r\nACGT\nAC\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineWidth(6, 5))
        ));

        let data = b">sq0\nACGT\nACGT\r\n";
        let mut indexer = Indexer::new(&data[..]);

        assert!(matches!(
            indexer.index_record(),
            Err(IndexError::InvalidLineWidth(6, 5))
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
