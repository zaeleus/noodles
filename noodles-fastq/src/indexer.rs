use std::{
    io::{self, BufRead},
    str,
};

use super::fai::Record;

/// A FASTQ indexer.
#[derive(Debug)]
pub struct Indexer<R> {
    inner: R,
    offset: u64,
    record: crate::Record,
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
            record: crate::Record::default(),
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
    /// let expected = fai::Record::new("r0", 4, 4, 4, 5, 11);
    /// assert_eq!(actual, Some(expected));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn index_record(&mut self) -> io::Result<Option<Record>> {
        use crate::io::reader::record::read_definition;

        // read name
        self.record.clear();
        self.offset += match read_definition(&mut self.inner, self.record.definition_mut()) {
            Ok(0) => return Ok(None),
            Ok(n) => n as u64,
            Err(e) => return Err(e),
        };

        let name = str::from_utf8(self.record.name())
            .map(String::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        // sequence
        let line_buf = self.record.name_mut();

        let sequence_offset = self.offset;

        line_buf.clear();
        self.offset += read_line(&mut self.inner, line_buf)? as u64;
        let line_width = line_buf.len() as u64;

        let line_bases = len_with_right_trim(line_buf) as u64;

        // plus line
        line_buf.clear();
        self.offset += read_line(&mut self.inner, line_buf)? as u64;

        // quality scores
        let quality_scores_offset = self.offset;

        line_buf.clear();
        self.offset += read_line(&mut self.inner, line_buf)? as u64;

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

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
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
@r1 LN:4
NNNNNNNNNN
+
NDLSNDLSND
";

        let mut indexer = Indexer::new(&data[..]);

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new("r0", 4, 4, 4, 5, 11)));

        let record = indexer.index_record()?;
        assert_eq!(record, Some(Record::new("r1", 10, 25, 10, 11, 38)));

        assert!(indexer.index_record()?.is_none());

        Ok(())
    }
}
