//! SAM reader and iterators.

pub mod record;
mod records;

pub use self::records::Records;

use std::io::{self, BufRead, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_fasta as fasta;

use self::record::Fields;
use super::{AlignmentReader, AlignmentRecord, Header};
use crate::lazy;

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
        self.records_with_fields(Fields::all())
    }

    /// Returns an iterator over records starting from the current stream position and decodes only
    /// the given fields.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::{self as sam, reader::record::Fields};
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut records = reader.records_with_fields(Fields::FLAGS);
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records_with_fields(&mut self, fields: Fields) -> Records<'_, R> {
        Records::new(self, fields)
    }

    /// Reads a single record without eagerly decoding its fields.
    ///
    /// This reads SAM fields from the underlying stream into the given record's buffer until a
    /// newline is reached. No fields are decoded, meaning the record is not necessarily valid.
    /// However, the structure of the byte stream is guaranteed to be record-like.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
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
    /// let mut record = sam::lazy::Record::default();
    /// reader.read_lazy_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        read_lazy_record(&mut self.inner, record)
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
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn AlignmentRecord>>> + 'a> {
        Box::new(
            self.records()
                .map(|result| result.map(|record| Box::new(record) as Box<dyn AlignmentRecord>)),
        )
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

fn read_lazy_record<R>(reader: &mut R, record: &mut lazy::Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.buf.clear();

    let mut len = 0;

    len += read_field(reader, &mut record.buf)?;
    record.bounds.read_name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.flags_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.reference_sequence_name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.alignment_start_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mapping_quality_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.cigar_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mate_reference_sequence_name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mate_alignment_start_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.template_length_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.sequence_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.quality_scores_end = record.buf.len();

    len += read_byte_line(reader, &mut record.buf)?;

    Ok(len)
}

fn read_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const DELIMITER: u8 = b'\t';

    let mut is_delimiter = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if is_delimiter || src.is_empty() {
            break;
        }

        let n = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                dst.extend_from_slice(&src[..i]);
                is_delimiter = true;
                i + 1
            }
            None => {
                dst.extend_from_slice(src);
                src.len()
            }
        };

        len += n;

        reader.consume(n);
    }

    Ok(len)
}

fn read_byte_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    match reader.read_until(LINE_FEED, buf)? {
        0 => Ok(0),
        n => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
                }
            }

            Ok(n)
        }
    }
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
