use std::io::{self, BufRead, Read};

use super::Record;

/// A FASTQ reader.
pub struct Reader<R> {
    inner: R,
    line_buf: Vec<u8>,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a FASTQ reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let reader = fastq::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            line_buf: Vec::new(),
        }
    }

    /// Reads a FASTQ record.
    ///
    /// This reads from the underlying stream until four lines are read: the read name, the
    /// sequence, the plus line, and the quality scores. Each line omits the trailing newline.
    ///
    /// The stream is expected to be at the start of a record.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::Reader::new(&data[..]);
    ///
    /// let mut record = fastq::Record::default();
    /// reader.read_record(&mut record)?;
    ///
    /// assert_eq!(record.read_name(), b"r0");
    /// assert_eq!(record.sequence(), b"ATCG");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        record.clear();

        let mut len = read_read_name(&mut self.inner, record.read_name_mut())?;

        if len > 0 {
            len += read_line(&mut self.inner, record.sequence_mut())?;

            self.line_buf.clear();
            len += read_line(&mut self.inner, &mut self.line_buf)?;

            len += read_line(&mut self.inner, record.quality_scores_mut())?;
        }

        Ok(len)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    let result = reader.read_until(b'\n', buf);
    buf.pop();
    result
}

fn read_read_name<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match consume_byte(reader, b'@') {
        Ok(n) => read_line(reader, buf).map(|m| m + n),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

fn consume_byte<R>(reader: &mut R, value: u8) -> io::Result<usize>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;

    if buf[0] == value {
        Ok(buf.len())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "read name missing @ prefix",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use crate::Record;

    use super::{read_line, Reader};

    #[test]
    fn test_read_record() {
        let data = "\
@noodles:1/1
AGCT
+
abcd
@noodles:2/1
TCGA
+
dcba
";

        let mut reader = Reader::new(data.as_bytes());
        let mut record = Record::default();

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 25);
        assert_eq!(record, Record::new("noodles:1/1", "AGCT", "abcd"));

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 25);
        assert_eq!(record, Record::new("noodles:2/1", "TCGA", "dcba"));

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 0);
    }

    #[test]
    fn test_read_line() {
        let data = "@fqlib\nAGCT\n";
        let mut reader = BufReader::new(data.as_bytes());

        let mut buf = Vec::new();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 7);
        assert_eq!(buf, b"@fqlib");

        buf.clear();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 5);
        assert_eq!(buf, b"AGCT");

        buf.clear();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 0);
    }
}
