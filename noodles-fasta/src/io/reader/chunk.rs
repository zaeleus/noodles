//! Chunk-based FASTA reader for high-throughput parsing.
//!
//! This reader processes FASTA data in large blocks, finding structural
//! delimiters (`>` and `\n`) in bulk using [`memchr`], and yields records
//! as zero-copy slices into the internal buffer.

use std::io::{self, Read};

use memchr::memchr;

const DEFAULT_BUF_SIZE: usize = 256 * 1024;

/// A borrowed FASTA record referencing data in the chunk reader's buffer.
///
/// The sequence may span multiple lines in the input. Use
/// [`sequence_lines()`](ChunkRecord::sequence_lines) for zero-copy access to
/// each line, or [`copy_sequence()`](ChunkRecord::copy_sequence) to
/// concatenate into a contiguous buffer.
#[derive(Clone, Debug)]
pub struct ChunkRecord<'a> {
    name: &'a [u8],
    description: &'a [u8],
    sequence_data: &'a [u8],
}

impl<'a> ChunkRecord<'a> {
    /// Returns the sequence name.
    pub fn name(&self) -> &'a [u8] {
        self.name
    }

    /// Returns the description.
    pub fn description(&self) -> &'a [u8] {
        self.description
    }

    /// Returns an iterator over the sequence lines (without newlines).
    pub fn sequence_lines(&self) -> impl Iterator<Item = &'a [u8]> {
        self.sequence_data
            .split(|&b| b == b'\n')
            .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
            .filter(|line| !line.is_empty())
    }

    /// Returns the total sequence length (excluding newlines).
    pub fn sequence_len(&self) -> usize {
        let whitespace = memchr::memchr2_iter(b'\n', b'\r', self.sequence_data).count();
        self.sequence_data.len() - whitespace
    }

    /// Copies the sequence into `buf`, stripping embedded newlines.
    /// Returns the number of bases copied.
    pub fn copy_sequence(&self, buf: &mut Vec<u8>) -> usize {
        let start = buf.len();
        for line in self.sequence_lines() {
            buf.extend_from_slice(line);
        }
        buf.len() - start
    }
}

/// A chunk-based FASTA reader that processes data in large blocks.
///
/// Instead of reading line-by-line, this reader fills a large internal buffer
/// and scans for record boundaries (`>`) and newlines in bulk using
/// SIMD-accelerated [`memchr`]. Records are returned as borrowed slices into
/// the buffer, avoiding per-record allocation.
///
/// # Examples
///
/// ```ignore
/// use noodles_fasta::io::reader::chunk::ChunkReader;
///
/// let data = b">seq0\nACGT\nTGCA\n>seq1\nAAAA\n";
/// let mut reader = ChunkReader::new(&data[..]);
///
/// let mut count = 0u64;
/// let mut total_bases = 0u64;
/// while let Some(record) = reader.next_record().unwrap() {
///     count += 1;
///     total_bases += record.sequence_len() as u64;
/// }
/// assert_eq!(count, 2);
/// assert_eq!(total_bases, 12);
/// ```
pub struct ChunkReader<R> {
    inner: R,
    buf: Vec<u8>,
    filled: usize,
    consumed: usize,
    finished: bool,
}

impl<R: Read> ChunkReader<R> {
    /// Creates a new chunk reader with the default buffer size (256 KiB).
    pub fn new(inner: R) -> Self {
        Self::with_capacity(DEFAULT_BUF_SIZE, inner)
    }

    /// Creates a new chunk reader with the specified buffer capacity.
    pub fn with_capacity(capacity: usize, inner: R) -> Self {
        Self {
            inner,
            buf: vec![0u8; capacity],
            filled: 0,
            consumed: 0,
            finished: false,
        }
    }

    /// Returns the next FASTA record, or `None` at EOF.
    pub fn next_record(&mut self) -> io::Result<Option<ChunkRecord<'_>>> {
        // Ensure we have data.
        if self.consumed >= self.filled {
            if self.finished {
                return Ok(None);
            }
            self.fill_buf()?;
            if self.filled == 0 {
                return Ok(None);
            }
        }

        // Skip any leading whitespace (between records).
        self.skip_whitespace();
        if self.consumed >= self.filled {
            if self.finished {
                return Ok(None);
            }
            self.compact_and_fill()?;
            self.skip_whitespace();
            if self.consumed >= self.filled {
                return Ok(None);
            }
        }

        loop {
            // Try to find a complete record in the buffer.
            if let Some(record_end) = self.find_record_end() {
                let record = self.parse_record(record_end);
                return Ok(Some(record));
            }

            if self.finished {
                // The rest of the buffer is the final record.
                let remaining = self.filled - self.consumed;
                if remaining == 0 {
                    return Ok(None);
                }
                let record_end = self.filled;
                let record = self.parse_record(record_end);
                return Ok(Some(record));
            }

            self.compact_and_fill()?;
        }
    }

    fn fill_buf(&mut self) -> io::Result<()> {
        self.consumed = 0;
        self.filled = 0;
        let n = read_full(&mut self.inner, &mut self.buf)?;
        self.filled = n;
        if n < self.buf.len() {
            self.finished = true;
        }
        Ok(())
    }

    fn compact_and_fill(&mut self) -> io::Result<()> {
        let remaining = self.filled - self.consumed;

        if remaining == 0 {
            self.consumed = 0;
            self.filled = 0;
        } else if self.consumed > 0 {
            self.buf.copy_within(self.consumed..self.filled, 0);
            self.consumed = 0;
            self.filled = remaining;
        } else {
            // Buffer is full but we can't find a record boundary — grow it.
            self.buf.resize(self.buf.len() * 2, 0);
        }

        let n = read_full(&mut self.inner, &mut self.buf[self.filled..])?;
        self.filled += n;
        if n == 0 || self.filled < self.buf.len() {
            self.finished = true;
        }
        Ok(())
    }

    fn skip_whitespace(&mut self) {
        while self.consumed < self.filled {
            match self.buf[self.consumed] {
                b'\n' | b'\r' | b' ' | b'\t' => self.consumed += 1,
                _ => break,
            }
        }
    }

    /// Find the end of the current record. A record starts with `>` and ends
    /// just before the next `>` that appears at the start of a line (i.e.,
    /// after a newline).
    fn find_record_end(&self) -> Option<usize> {
        let data = &self.buf[self.consumed..self.filled];

        // The current position should be at '>'. Skip past the header line first.
        let header_end = memchr(b'\n', data)?;
        let mut pos = header_end + 1;

        // Scan for '>' directly using memchr, then verify it's at line start.
        while pos < data.len() {
            match memchr(b'>', &data[pos..]) {
                Some(i) => {
                    let abs = pos + i;
                    // Valid record boundary: '>' preceded by '\n'.
                    if abs > 0 && data[abs - 1] == b'\n' {
                        return Some(self.consumed + abs);
                    }
                    // Not at line start (e.g., '>' in a header description).
                    // Skip past this '>' and keep scanning.
                    pos = abs + 1;
                }
                None => return None,
            }
        }

        None
    }

    /// Parse a record from `self.consumed` to `record_end`.
    fn parse_record(&mut self, record_end: usize) -> ChunkRecord<'_> {
        let data = &self.buf[self.consumed..record_end];
        self.consumed = record_end;

        // Find the header line end.
        let header_end = memchr(b'\n', data).unwrap_or(data.len());

        // Parse definition: >name [description]
        let definition_line = &data[1..header_end]; // skip '>'
        let definition_line = definition_line
            .strip_suffix(b"\r")
            .unwrap_or(definition_line);
        let (name, description) = split_definition(definition_line);

        // Everything after the header line is sequence data (with embedded newlines).
        let seq_start = if header_end < data.len() {
            header_end + 1
        } else {
            data.len()
        };
        let sequence_data = &data[seq_start..];
        // Strip trailing newlines from the sequence block.
        let sequence_data = sequence_data
            .strip_suffix(b"\r\n")
            .or_else(|| sequence_data.strip_suffix(b"\n"))
            .unwrap_or(sequence_data);

        ChunkRecord {
            name,
            description,
            sequence_data,
        }
    }
}

fn split_definition(line: &[u8]) -> (&[u8], &[u8]) {
    if let Some(i) = memchr::memchr2(b' ', b'\t', line) {
        (&line[..i], &line[i + 1..])
    } else {
        (line, &[])
    }
}

fn read_full(reader: &mut impl Read, mut buf: &mut [u8]) -> io::Result<usize> {
    let mut total = 0;
    while !buf.is_empty() {
        match reader.read(buf) {
            Ok(0) => break,
            Ok(n) => {
                total += n;
                buf = &mut buf[n..];
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_records() {
        let data = b">seq0\nACGT\n>seq1\nTGCA\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq0");
        assert_eq!(record.description(), b"");
        assert_eq!(record.sequence_len(), 4);
        let mut seq = Vec::new();
        record.copy_sequence(&mut seq);
        assert_eq!(seq, b"ACGT");

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq1");
        let mut seq = Vec::new();
        record.copy_sequence(&mut seq);
        assert_eq!(seq, b"TGCA");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_multiline_sequence() {
        let data = b">seq0\nACGT\nTGCA\nAAAA\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.sequence_len(), 12);
        let mut seq = Vec::new();
        record.copy_sequence(&mut seq);
        assert_eq!(seq, b"ACGTTGCAAAAA");
    }

    #[test]
    fn test_with_description() {
        let data = b">seq0 some description\nACGT\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq0");
        assert_eq!(record.description(), b"some description");
    }

    #[test]
    fn test_empty_input() {
        let data = b"";
        let mut reader = ChunkReader::new(&data[..]);
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_boundary_crossing() {
        let data = b">seq0\nACGTACGTACGT\n>seq1\nTGCATGCATGCA\n";
        let mut reader = ChunkReader::with_capacity(16, &data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq0");
        assert_eq!(record.sequence_len(), 12);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq1");
        assert_eq!(record.sequence_len(), 12);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_crlf_line_endings() {
        let data = b">seq0\r\nACGT\r\nTGCA\r\n>seq1\r\nAAAA\r\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq0");
        assert_eq!(record.sequence_len(), 8);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq1");
        assert_eq!(record.sequence_len(), 4);
    }

    #[test]
    fn test_many_records_small_buffer() {
        let mut data = Vec::new();
        for i in 0..100 {
            data.extend_from_slice(format!(">seq_{i}\n").as_bytes());
            data.extend_from_slice(b"ACGTACGT\nTGCATGCA\n");
        }

        let mut reader = ChunkReader::with_capacity(64, &data[..]);
        let mut count = 0;
        let mut total_bases = 0;
        while let Some(record) = reader.next_record().unwrap() {
            count += 1;
            total_bases += record.sequence_len();
        }
        assert_eq!(count, 100);
        assert_eq!(total_bases, 1600);
    }

    #[test]
    fn test_no_trailing_newline() {
        let data = b">seq0\nACGT";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"seq0");
        assert_eq!(record.sequence_len(), 4);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_single_long_sequence() {
        let mut data = b">chr1\n".to_vec();
        // 10 lines of 80 bases each = 800 bases
        for _ in 0..10 {
            data.extend_from_slice(&[b'A'; 80]);
            data.push(b'\n');
        }

        let mut reader = ChunkReader::with_capacity(128, &data[..]);
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"chr1");
        assert_eq!(record.sequence_len(), 800);
    }
}
