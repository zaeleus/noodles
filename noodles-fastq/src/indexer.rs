use std::io::{self, BufRead};

use super::fai::Record;

const LINE_FEED: u8 = b'\n';

/// A FASTQ indexer.
#[derive(Debug)]
pub struct Indexer<R> {
    inner: R,
    offset: u64,
    line_buf: Vec<u8>,
}

impl<R> Indexer<R>
where
    R: BufRead,
{
    /// Creates a FASTQ indexer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = b"@r0\nACTG\n+\nNDLS\n";
    /// let mut indexer = fastq::Indexer::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            offset: 0,
            line_buf: Vec::new(),
        }
    }

    /// Indexes a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq::{self as fastq, fai};
    ///
    /// let data = b"@r0\nACTG\n+\nNDLS\n";
    /// let mut indexer = fastq::Indexer::new(&data[..]);
    ///
    /// let actual = indexer.index_record()?;
    /// let expected = fai::Record::new(String::from("r0"), 4, 4, 4, 5, 11);
    /// assert_eq!(actual, Some(expected));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn index_record(&mut self) -> io::Result<Option<Record>> {
        // read name
        self.line_buf.clear();
        self.offset += match read_name(&mut self.inner, &mut self.line_buf) {
            Ok(0) => return Ok(None),
            Ok(n) => n as u64,
            Err(e) => return Err(e),
        };

        let name = String::from_utf8(self.line_buf.clone())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        // sequence
        let sequence_offset = self.offset;

        self.line_buf.clear();
        self.offset += read_line(&mut self.inner, &mut self.line_buf)? as u64;
        let line_width = self.line_buf.len() as u64;

        let line_bases = len_with_right_trim(&self.line_buf) as u64;

        // plus line
        self.line_buf.clear();
        self.offset += read_line(&mut self.inner, &mut self.line_buf)? as u64;

        // quality scores
        let quality_scores_offset = self.offset;

        self.line_buf.clear();
        self.offset += read_line(&mut self.inner, &mut self.line_buf)? as u64;

        Ok(Some(Record::new(
            name,
            line_bases,
            sequence_offset,
            line_bases,
            line_width,
            quality_scores_offset,
        )))
    }
}

fn read_name<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    use crate::reader::read_u8;

    const NAME_PREFIX: u8 = b'@';

    match read_u8(reader) {
        Ok(NAME_PREFIX) => crate::reader::read_line(reader, buf).map(|n| n + 1),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid name prefix",
        )),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    reader.read_until(LINE_FEED, buf)
}

fn len_with_right_trim(buf: &[u8]) -> usize {
    match buf.iter().rposition(|b| !b.is_ascii_whitespace()) {
        Some(i) => i + 1,
        None => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_record() -> io::Result<()> {
        let data = b"\
@r0
ACGT
+
NDLS
@r1
NNNNNNNNNN
+
NDLSNDLSND
";

        let mut indexer = Indexer::new(&data[..]);

        let record = indexer.index_record()?;
        assert_eq!(
            record,
            Some(Record::new(String::from("r0"), 4, 4, 4, 5, 11))
        );

        let record = indexer.index_record()?;
        assert_eq!(
            record,
            Some(Record::new(String::from("r1"), 10, 20, 10, 11, 33))
        );

        assert!(indexer.index_record()?.is_none());

        Ok(())
    }
}
