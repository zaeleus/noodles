//! VCF reader and iterators.

mod query;
mod records;

pub use self::{query::Query, records::Records};

use std::io::{self, BufRead, Read, Seek};

use memchr::memchr;
use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Region};
use noodles_csi::BinningIndex;
use noodles_tabix as tabix;

use super::Header;

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

const HEADER_PREFIX: u8 = b'#';

/// A VCF reader.
///
/// The VCF format has two main parts: 1) a header and 2) a list of VCF records.
///
/// Each header line is prefixed with a `#` (number sign) and is terminated by the header header
/// (`#CHROM`...; inclusive).
///
/// VCF records are line-based and follow directly after the header until EOF.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io::BufReader};
/// use noodles_vcf as vcf;
///
/// let mut reader = File::open("sample.vcf")
///     .map(BufReader::new)
///     .map(vcf::Reader::new)?;
///
/// let header = reader.read_header()?.parse()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     println!("{:?}", record);
/// }
/// # Ok::<_, Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a VCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let reader = vcf::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let reader = vcf::Reader::new(&data[..]);
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
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let reader = vcf::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the raw VCF header.
    ///
    /// This reads all header lines prefixed with a `#` (number sign), which includes the header
    /// header (`#CHROM`...).
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw VCF header as a [`String`], and as such, it is not necessarily valid.
    /// The raw header can subsequently be parsed as a [`crate::Header`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// assert_eq!(header, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        read_header(&mut self.inner)
    }

    /// Reads a single raw VCF record.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline character. The buffer can subsequently be parsed as a
    /// [`crate::Record`].
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
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut buf = String::new();
    /// reader.read_record(&mut buf)?;
    ///
    /// assert_eq!(buf, "sq0\t1\t.\tA\t.\t.\tPASS\t.");
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
    /// ```no_run
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// let header = reader.read_header()?.parse()?;
    ///
    /// let mut records = reader.records(&header);
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn records<'r, 'h>(&'r mut self, header: &'h Header) -> Records<'r, 'h, R> {
        Records::new(self, header)
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let data = Vec::new();
    /// let reader = vcf::Reader::new(bgzf::Reader::new(&data[..]));
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(virtual_position.compressed(), 0);
    /// assert_eq!(virtual_position.uncompressed(), 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
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
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let data = Cursor::new(Vec::new());
    /// let mut reader = vcf::Reader::new(bgzf::Reader::new(data));
    ///
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }

    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bgzf as bgzf;;
    /// use noodles_core::Region;
    /// use noodles_tabix as tabix;
    /// use noodles_vcf as vcf;
    ///
    /// let mut reader = File::open("sample.vcf.gz")
    ///     .map(bgzf::Reader::new)
    ///     .map(vcf::Reader::new)?;
    ///
    /// let header = reader.read_header()?.parse()?;
    ///
    /// let index = tabix::read("sample.vcf.gz.tbi")?;
    /// let region = Region::mapped("sq0", 8..=13);
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h>(
        &'r mut self,
        header: &'h Header,
        index: &tabix::Index,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, R>> {
        let (reference_sequence_id, reference_sequence_name, interval) =
            resolve_region(index, region)?;

        let chunks = index.query(reference_sequence_id, interval)?;

        Ok(Query::new(
            self,
            chunks,
            reference_sequence_name,
            interval,
            header,
        ))
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

        let (read_eol, len) = if let Some(i) = memchr(LINE_FEED as u8, buf) {
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

// Reads all bytes until a line feed ('\n') or EOF is reached.
//
// The buffer will not include the trailing newline ('\n' or '\r\n').
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

pub(crate) fn resolve_region(
    index: &tabix::Index,
    region: &Region,
) -> io::Result<(usize, String, Interval)> {
    if let Some(r) = region.as_mapped() {
        let i = index
            .reference_sequence_names()
            .get_index_of(r.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "region reference sequence does not exist in reference sequences: {:?}",
                        region
                    ),
                )
            })?;

        Ok((i, r.name().into(), r.interval()))
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "region is not mapped",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;

    static DATA: &[u8] = b"\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t8
sq0\t13
";

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = DATA;

        let actual = read_header(&mut reader)?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_no_records() -> io::Result<()> {
        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let mut reader = expected.as_bytes();
        let actual = read_header(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let mut reader = BufReader::with_capacity(16, expected.as_bytes());
        let actual = read_header(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_no_header() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;
        assert!(actual.is_empty());

        let data = b"sq0\t1\t.\tA\t.\t.\tPASS\t.\n";
        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;
        assert!(actual.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut reader = Reader::new(DATA);
        reader.read_header()?;

        let mut buf = String::new();
        let bytes_read = reader.read_record(&mut buf)?;
        assert_eq!(bytes_read, 6);
        assert_eq!(buf, "sq0\t8");

        buf.clear();
        let bytes_read = reader.read_record(&mut buf)?;
        assert_eq!(bytes_read, 7);
        assert_eq!(buf, "sq0\t13");

        buf.clear();
        let bytes_read = reader.read_record(&mut buf)?;
        assert_eq!(bytes_read, 0);
        assert!(buf.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_line() -> io::Result<()> {
        let mut buf = String::new();

        let data = b"noodles\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        let data = b"noodles\r\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        let data = b"noodles";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        Ok(())
    }
}
