//! Chunk-based FASTQ reader for high-throughput parsing.
//!
//! This reader processes FASTQ data in large blocks, finding all newline
//! positions in bulk using [`memchr`], and yields records as zero-copy slices
//! into the internal buffer.

use std::io::{self, Read};

use memchr::memchr;

const DEFAULT_BUF_SIZE: usize = 256 * 1024;

/// A borrowed FASTQ record referencing data in the chunk reader's buffer.
#[derive(Clone, Debug)]
pub struct ChunkRecord<'a> {
    name: &'a [u8],
    description: &'a [u8],
    sequence: &'a [u8],
    quality_scores: &'a [u8],
}

impl<'a> ChunkRecord<'a> {
    /// Returns the read name.
    pub fn name(&self) -> &'a [u8] {
        self.name
    }

    /// Returns the description.
    pub fn description(&self) -> &'a [u8] {
        self.description
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    /// Returns the quality scores.
    pub fn quality_scores(&self) -> &'a [u8] {
        self.quality_scores
    }
}

/// A chunk-based FASTQ reader that processes data in large blocks.
///
/// Instead of reading line-by-line, this reader fills a large internal buffer
/// and scans for all newline positions at once using SIMD-accelerated
/// [`memchr`]. Records are returned as borrowed slices into the buffer,
/// avoiding per-record allocation.
///
/// # Examples
///
/// ```ignore
/// use noodles_fastq::io::reader::chunk::ChunkReader;
///
/// let data = b"@r0\nACGT\n+\nNDLS\n@r1\nTGCA\n+\nSLDN\n";
/// let mut reader = ChunkReader::new(&data[..]);
///
/// let mut count = 0u64;
/// let mut total_bases = 0u64;
/// while let Some(record) = reader.next_record().unwrap() {
///     count += 1;
///     total_bases += record.sequence().len() as u64;
/// }
/// assert_eq!(count, 2);
/// assert_eq!(total_bases, 8);
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

    /// Returns the next FASTQ record, or `None` at EOF.
    pub fn next_record(&mut self) -> io::Result<Option<ChunkRecord<'_>>> {
        // Ensure we have data in the buffer.
        if self.consumed >= self.filled {
            if self.finished {
                return Ok(None);
            }
            self.fill_buf()?;
            if self.filled == 0 {
                return Ok(None);
            }
        }

        loop {
            // Try to find a complete record in the current buffer.
            if let Some(record_end) = self.find_record_end() {
                let record = self.parse_record(record_end);
                return Ok(Some(record));
            }

            // Not enough data for a complete record. If we're at EOF, that's
            // an error (truncated record) unless we have no data at all.
            if self.finished {
                let remaining = self.filled - self.consumed;
                if remaining == 0 {
                    return Ok(None);
                }
                // Check if remaining data is just whitespace.
                let leftover = &self.buf[self.consumed..self.filled];
                if leftover.iter().all(|&b| b == b'\n' || b == b'\r') {
                    return Ok(None);
                }
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated FASTQ record",
                ));
            }

            // Compact the buffer and read more data.
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
            // Buffer is full but we can't find a record — grow it.
            self.buf.resize(self.buf.len() * 2, 0);
        }

        let n = read_full(&mut self.inner, &mut self.buf[self.filled..])?;
        self.filled += n;
        if n == 0 || self.filled < self.buf.len() {
            self.finished = true;
        }
        Ok(())
    }

    /// Find the end position (exclusive) of the next complete record.
    /// A FASTQ record has exactly 4 lines, so we need to find 4 newlines.
    fn find_record_end(&self) -> Option<usize> {
        let data = &self.buf[self.consumed..self.filled];
        let mut pos = 0;

        for _ in 0..4 {
            let remaining = &data[pos..];
            match memchr(b'\n', remaining) {
                Some(i) => pos += i + 1,
                None => return None,
            }
        }

        Some(self.consumed + pos)
    }

    /// Parse a record from `self.consumed` to `record_end`.
    /// Updates `self.consumed` to `record_end`.
    fn parse_record(&mut self, record_end: usize) -> ChunkRecord<'_> {
        let data = &self.buf[self.consumed..record_end];
        self.consumed = record_end;

        // Find the 4 newline positions within this record's data.
        let mut newlines = [0usize; 4];
        let mut pos = 0;
        for nl in &mut newlines {
            *nl = pos + memchr(b'\n', &data[pos..]).unwrap();
            pos = *nl + 1;
        }

        // Line 1: @name [description]\n
        let definition_line = &data[1..newlines[0]]; // skip '@'
        let definition_line = strip_cr(definition_line);
        let (name, description) = split_definition(definition_line);

        // Line 2: sequence\n
        let sequence = strip_cr(&data[newlines[0] + 1..newlines[1]]);

        // Line 3: +[ignored]\n (skip)
        // Line 4: quality\n
        let quality_scores = strip_cr(&data[newlines[2] + 1..newlines[3]]);

        ChunkRecord {
            name,
            description,
            sequence,
            quality_scores,
        }
    }
}

fn strip_cr(s: &[u8]) -> &[u8] {
    s.strip_suffix(b"\r").unwrap_or(s)
}

fn split_definition(line: &[u8]) -> (&[u8], &[u8]) {
    if let Some(i) = memchr::memchr2(b' ', b'\t', line) {
        (&line[..i], &line[i + 1..])
    } else {
        (line, &[])
    }
}

/// Read as much as possible into `buf`, retrying on short reads.
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
        let data = b"@r0\nACGT\n+\nNDLS\n@r1\nTGCA\n+\nSLDN\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.description(), b"");
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r1");
        assert_eq!(record.sequence(), b"TGCA");
        assert_eq!(record.quality_scores(), b"SLDN");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_with_description() {
        let data = b"@r0 some description\nACGT\n+\nNDLS\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.description(), b"some description");
        assert_eq!(record.sequence(), b"ACGT");
    }

    #[test]
    fn test_empty_input() {
        let data = b"";
        let mut reader = ChunkReader::new(&data[..]);
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_boundary_crossing() {
        // Use a tiny buffer to force records to span boundaries.
        let data = b"@r0\nACGTACGT\n+\nNDLSNDLS\n@r1\nTGCATGCA\n+\nSLDNSLDN\n";
        let mut reader = ChunkReader::with_capacity(16, &data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.sequence(), b"ACGTACGT");
        assert_eq!(record.quality_scores(), b"NDLSNDLS");

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r1");
        assert_eq!(record.sequence(), b"TGCATGCA");
        assert_eq!(record.quality_scores(), b"SLDNSLDN");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_crlf_line_endings() {
        let data = b"@r0\r\nACGT\r\n+\r\nNDLS\r\n";
        let mut reader = ChunkReader::new(&data[..]);

        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");
    }

    #[test]
    fn test_many_records_small_buffer() {
        let mut data = Vec::new();
        for i in 0..100 {
            data.extend_from_slice(format!("@read_{i}\n").as_bytes());
            data.extend_from_slice(b"ACGTACGT\n+\nFFFFFFFF\n");
        }

        let mut reader = ChunkReader::with_capacity(64, &data[..]);
        let mut count = 0;
        while reader.next_record().unwrap().is_some() {
            count += 1;
        }
        assert_eq!(count, 100);
    }

    #[test]
    fn test_truncated_record() {
        let data = b"@r0\nACGT\n+\n";
        let mut reader = ChunkReader::new(&data[..]);
        assert!(reader.next_record().is_err());
    }
}
