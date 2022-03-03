//! SAM reader and iterators.

mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_fasta as fasta;

use super::{AlignmentReader, Header, Record};

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

const HEADER_PREFIX: u8 = b'@';

/// A SAM reader.
///
/// The SAM format is comprised of two parts: 1) a header and 2) a list of records.
///
/// Each header line is prefixed with an `@` (at sign). The header is optional and may be empty.
///
/// SAM records are line-based and follow directly after the header or the start of the file until
/// EOF.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io::{self, BufReader}};
/// use noodles_sam as sam;
///
/// let mut reader = File::open("sample.sam")
///     .map(BufReader::new)
///     .map(sam::Reader::new)?;
///
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
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let reader = sam::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::Reader::new(&data[..]);
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
    /// use noodles_sam as sam;
    /// let data = [];
    /// let mut reader = sam::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the raw SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`crate::Header`].
    ///
    /// The SAM header is optional, and if it is missing, an empty string is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// assert_eq!(header, "@HD\tVN:1.6\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        read_header(&mut self.inner)
    }

    /// Reads a single raw SAM record.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline character.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`]), but using
    /// this method allows control of the line buffer and whether the raw record should be parsed.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut buf = String::new();
    /// reader.read_record(&mut buf)?;
    /// assert_eq!(buf, "*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*");
    ///
    /// assert_eq!(reader.read_record(&mut buf)?, 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// Unlike [`Self::read_record`], each record is parsed as a [`crate::Record`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut records = reader.records();
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read + Seek,
{
    /// Seeks the underlying BGZF stream to the given virtual position.
    ///
    /// Virtual positions typically come from an associated index.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use std::io::BufReader;
    ///
    /// use noodles_bgzf as bgzf;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.sam.gz")
    ///     .map(bgzf::Reader::new)
    ///     .map(sam::Reader::new)?;
    ///
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }
}

impl<R> AlignmentReader for Reader<R>
where
    R: BufRead,
{
    fn read_alignment_header(&mut self) -> io::Result<Header> {
        read_header(&mut self.inner).and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    fn alignment_records<'a>(
        &'a mut self,
        _: &'a fasta::Repository,
        _: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<Record>> + 'a> {
        Box::new(self.records())
    }
}

fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: BufRead,
{
    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = reader.fill_buf()?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != HEADER_PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = buf.iter().position(|&b| b == LINE_FEED as u8) {
            header_buf.extend(&buf[..=i]);
            (true, i + 1)
        } else {
            header_buf.extend(buf);
            (false, buf.len())
        };

        is_eol = read_eol;

        reader.consume(len);
    }

    String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_line(buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;

    #[test]
    fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();
        let header = read_header(&mut reader)?;
        assert_eq!(header, data);
        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());
        let header = read_header(&mut reader)?;
        assert_eq!(header, data);
        Ok(())
    }

    #[test]
    fn test_read_line() -> io::Result<()> {
        fn t(buf: &mut String, mut reader: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = String::new();

        t(&mut buf, b"noodles\n", "noodles")?;
        t(&mut buf, b"noodles\r\n", "noodles")?;
        t(&mut buf, b"noodles", "noodles")?;

        Ok(())
    }
}
