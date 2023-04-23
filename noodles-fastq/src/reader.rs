//! FASTQ reader.

mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Read};

use super::Record;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

/// A FASTQ reader.
pub struct Reader<R> {
    inner: R,
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
        Self { inner }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = [];
    /// let reader = fastq::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = [];
    /// let mut reader = fastq::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = [];
    /// let reader = fastq::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
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
    /// assert_eq!(record.name(), b"r0");
    /// assert_eq!(record.sequence(), b"ATCG");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record)
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

fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.clear();

    let mut len = match read_definition(reader, record) {
        Ok(0) => return Ok(0),
        Ok(n) => n,
        Err(e) => return Err(e),
    };

    len += read_line(reader, record.sequence_mut())?;
    len += consume_description(reader)?;
    len += read_line(reader, record.quality_scores_mut())?;

    Ok(len)
}

pub(crate) fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
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

fn consume_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr;

    let mut is_eol = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if src.is_empty() || is_eol {
            break;
        }

        let n = match memchr(LINE_FEED, src) {
            Some(i) => {
                is_eol = true;
                i + 1
            }
            None => src.len(),
        };

        reader.consume(n);

        len += n;
    }

    Ok(len)
}

pub(super) fn read_definition<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr2;

    const DELIMITER: u8 = b' ';
    const NAME_PREFIX: u8 = b'@';

    match read_u8(reader) {
        Ok(prefix) => {
            if prefix != NAME_PREFIX {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid name prefix",
                ));
            }
        }
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    }

    let mut is_eol = false;
    let mut len = 1;

    loop {
        let src = reader.fill_buf()?;

        if src.is_empty() {
            break;
        }

        let (matched_needle, n) = match memchr2(DELIMITER, LINE_FEED, src) {
            Some(i) => {
                let name_src = match src[i] {
                    DELIMITER => &src[..i],
                    LINE_FEED => {
                        is_eol = true;

                        if src.ends_with(&[CARRIAGE_RETURN]) {
                            &src[..i - 1]
                        } else {
                            &src[..i]
                        }
                    }
                    _ => unreachable!(),
                };

                record.name_mut().extend(name_src);

                (true, i + 1)
            }
            None => {
                record.name_mut().extend(src);
                (false, src.len())
            }
        };

        len += n;
        reader.consume(n);

        if matched_needle {
            break;
        }
    }

    if !is_eol {
        len += read_line(reader, record.description_mut())?;
    }

    Ok(len)
}

fn consume_description<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    const DESCRIPTION_PREFIX: u8 = b'+';

    match read_u8(reader)? {
        DESCRIPTION_PREFIX => consume_line(reader).map(|n| n + 1),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid description prefix",
        )),
    }
}

fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles:2/1
TCGA
+noodles:2/1
dcba
";

        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        let expected = Record::new("noodles:1/1", "AGCT", "abcd");
        assert_eq!(record, expected);

        read_record(&mut reader, &mut record)?;
        let expected = Record::new("noodles:2/1", "TCGA", "dcba");
        assert_eq!(record, expected);

        let n = read_record(&mut reader, &mut record)?;
        assert_eq!(n, 0);

        Ok(())
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

    #[test]
    fn test_read_definition() -> io::Result<()> {
        let mut record = Record::default();

        let data = b"@r0\n";
        let mut reader = &data[..];
        record.clear();
        read_definition(&mut reader, &mut record)?;
        assert_eq!(record.name(), b"r0");
        assert!(record.description().is_empty());

        let data = b"@r0 LN:4\n";
        let mut reader = &data[..];
        record.clear();
        read_definition(&mut reader, &mut record)?;
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.description(), b"LN:4");

        let data = b"r0\n";
        let mut reader = &data[..];
        record.clear();
        assert!(matches!(
            read_definition(&mut reader, &mut record),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_consume_description() -> io::Result<()> {
        let data = b"+r0\n";
        let mut reader = &data[..];
        consume_description(&mut reader)?;

        let data = b"r0\n";
        let mut reader = &data[..];
        assert!(matches!(
            consume_description(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
