mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Read};

use super::Record;

const READ_NAME_PREFIX: u8 = b'@';
const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

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

        let mut len = match read_read_name(&mut self.inner, record.read_name_mut()) {
            Ok(0) => return Ok(0),
            Ok(n) => n,
            Err(e) => return Err(e),
        };

        len += read_line(&mut self.inner, record.sequence_mut())?;

        self.line_buf.clear();
        len += read_line(&mut self.inner, &mut self.line_buf)?;

        len += read_line(&mut self.inner, record.quality_scores_mut())?;

        Ok(len)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a record.
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
    /// let mut records = reader.records();
    ///
    /// assert_eq!(
    ///     records.next().transpose()?,
    ///     Some(fastq::Record::new("r0", "ATCG", "NDLS")
    /// ));
    ///
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_until(LINE_FEED, buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

pub(crate) fn read_read_name<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match consume_byte(reader, READ_NAME_PREFIX) {
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
    use super::*;

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
    fn test_read_line() -> io::Result<()> {
        let mut buf = Vec::new();

        let data = b"noodles\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles\r\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        Ok(())
    }
}
