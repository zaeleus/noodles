//! SAM reader and iterators

mod records;

pub use self::records::Records;

use std::io::{self, BufRead};

const HEADER_PREFIX: u8 = b'@';
const NEWLINE: u8 = b'\n';

/// A SAM reader.
///
/// The SAM format is comprised to two parts:
///
///   1. a SAM header, and
///   2. a list of SAM records.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io::{self, BufReader}};
/// use noodles_sam as sam;
///
/// let mut reader = File::open("sample.sam").map(BufReader::new).map(sam::Reader::new)?;
/// reader.read_header()?;
///
/// for result in reader.records() {
///     let record = result?;
///     println!("{:?}", record);
/// }
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a SAM reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::{self, BufReader}};
    /// use noodles_sam as sam;
    /// let mut reader = File::open("sample.sam").map(BufReader::new).map(sam::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the raw SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`sam::Header`].
    ///
    /// [`String`]: https://doc.rust-lang.org/std/string/struct.String.html
    /// [`sam::Header`]: header/struct.Header.html
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::{self, BufReader}};
    /// use noodles_sam as sam;
    /// let mut reader = File::open("sample.sam").map(BufReader::new).map(sam::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        let mut header_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = self.inner.fill_buf()?;

            if eol && !buf.is_empty() && buf[0] != HEADER_PREFIX {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == NEWLINE) {
                Some(i) => {
                    header_buf.extend(&buf[..=i]);
                    (true, i + 1)
                }
                None => {
                    header_buf.extend(buf);
                    (false, buf.len())
                }
            };

            eol = read_eol;
            self.inner.consume(len);
        }

        String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Reads a single raw SAM record.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline character.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`records`]), but using this
    /// method allows control of the line buffer and whether the raw record should be parsed.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// [`records`]: #method.records
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::{self, BufReader}};
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.sam").map(BufReader::new).map(sam::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut record = String::new();
    /// reader.read_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        let result = self.inner.read_line(buf);
        buf.pop();
        result
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// Unlike [`read_record`], each record is parsed as a [`sam::Record`].
    ///
    /// [`read_record`]: #method.read_record
    /// [`sam::Record`]: record/struct.Record.html
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::{self, BufReader}};
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.sam").map(BufReader::new).map(sam::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static DATA: &[u8] = b"\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
r001
r002
";

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = Reader::new(&DATA[..]);

        let actual = reader.read_header()?;
        let expected = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
";
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut reader = Reader::new(&DATA[..]);
        reader.read_header()?;

        let mut buf = String::new();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "r001");

        buf.clear();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "r002");

        buf.clear();
        let len = reader.read_record(&mut buf)?;
        assert_eq!(len, 0);

        Ok(())
    }
}
