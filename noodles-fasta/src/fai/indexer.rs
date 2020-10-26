//! FASTA indexer.

use std::{
    error::Error,
    fmt,
    io::{self, BufRead},
};

use memchr::memchr;

use crate::{
    reader::{DEFINITION_PREFIX, NEWLINE},
    record::definition::{Definition, ParseError},
};

use super::Record;

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

    /// Reads a single sequence line.
    ///
    /// Trailing whitespaces are not discarded.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the entire sequence of the current record was read.
    fn read_sequence_line(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut bytes_read = 0;
        let mut is_eol = false;

        loop {
            let reader_buf = self.inner.fill_buf()?;

            if is_eol || reader_buf.is_empty() || reader_buf[0] == DEFINITION_PREFIX {
                break;
            }

            let len = match memchr(NEWLINE, reader_buf) {
                Some(i) => {
                    // Do not consume newline.
                    buf.extend(&reader_buf[..=i]);
                    is_eol = true;
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

    /// Indexes a raw FASTA record.
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// definition.
    ///
    /// # Errors
    ///
    /// An error is returned if the record fails to be completely read. This includes when
    ///
    ///   * the stream is not at the start of a definition,
    ///   * the record is missing a sequence,
    ///   * or the sequence lines are not the same length, excluding the last line.
    pub fn index_record(&mut self) -> Result<Option<Record>, IndexError> {
        let mut def_str = String::new();
        let mut line = Vec::new();
        let mut length = 0u64;

        let def_len = self.inner.read_line(&mut def_str)?;
        if def_len == 0 {
            return Ok(None);
        }

        self.offset += def_len as u64;
        let index_offset = self.offset;
        let def: Definition = def_str.trim_end().parse()?;

        // The first sequence line determines how long each line should be.
        let mut seq_len = self.read_sequence_line(&mut line)? as u64;
        let line_width = seq_len;
        let line_bases = len_with_right_trim(&line) as u64;

        loop {
            let prev_line_width = seq_len;
            let prev_line_bases = len_with_right_trim(&line) as u64;
            length += prev_line_bases;
            self.offset += prev_line_width;

            line.clear();
            seq_len = self.read_sequence_line(&mut line)? as u64;
            if seq_len == 0 {
                break;
            // If there are more lines, check the previous line has equal length to first.
            } else if prev_line_width != line_width || prev_line_bases != line_bases {
                return Err(IndexError::InvalidLineLength(self.offset));
            }
        }

        if length == 0 {
            return Err(IndexError::EmptySequence(self.offset));
        }

        let record = Record::new(
            def.reference_sequence_name().to_string(),
            length,
            index_offset,
            line_bases,
            line_width,
        );

        Ok(Some(record))
    }
}

fn len_with_right_trim(vec: &[u8]) -> usize {
    match vec.iter().rposition(|x| !x.is_ascii_whitespace()) {
        Some(i) => i + 1,
        None => 0,
    }
}

#[derive(Debug)]
pub enum IndexError {
    EmptySequence(u64),
    InvalidDefinition(ParseError),
    InvalidLineLength(u64),
    IoError(io::Error),
}

impl Error for IndexError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::EmptySequence(_) => None,
            Self::InvalidDefinition(e) => Some(e),
            Self::InvalidLineLength(_) => None,
            Self::IoError(e) => Some(e),
        }
    }
}

impl fmt::Display for IndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequence(offset) => write!(f, "empty sequence at offset {}", offset),
            Self::InvalidDefinition(e) => e.fmt(f),
            Self::InvalidLineLength(offset) => {
                write!(f, "different line lengths at offset {}", offset)
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
            _ => io::Error::new(io::ErrorKind::InvalidInput, error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_sequence_line() -> io::Result<()> {
        let data = b"ACGT\nNNNN\n";
        let mut indexer = Indexer::new(&data[..]);

        let mut buf = Vec::new();
        indexer.read_sequence_line(&mut buf)?;
        assert_eq!(buf, b"ACGT\n");

        Ok(())
    }

    #[test]
    fn test_index_record() -> Result<(), IndexError> {
        let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
        let mut indexer = Indexer::new(&data[..]);

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new(String::from("sq0"), 4, 5, 4, 5)));

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new(String::from("sq1"), 10, 15, 4, 5)));

        assert!(indexer.index_record()?.is_none());

        Ok(())
    }

    #[test]
    fn test_len_with_right_trim() {
        assert_eq!(len_with_right_trim(b"ATGC\n"), 4);
        assert_eq!(len_with_right_trim(b"ATGC\r\n"), 4);
    }
}
