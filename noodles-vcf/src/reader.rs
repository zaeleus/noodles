//! VCF reader and iterators.

mod records;

pub use self::records::Records;

use std::io::{self, BufRead, BufReader, Read, Seek};

use noodles_bgzf as bgzf;

const LINE_FEED: u8 = b'\n';
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
/// # use std::{fs::File, io::{self, BufReader}};
/// use noodles_vcf as vcf;
///
/// let mut reader = File::open("sample.vcf")
///     .map(BufReader::new)
///     .map(vcf::Reader::new)?;
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

    /// Reads the raw VCF header.
    ///
    /// This reads all header lines prefixed with a `#` (number sign) and is terminated by the
    /// header header (`#CHROM`...; inclusive).
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`vcf::Header`].
    ///
    /// [`String`]: https://doc.rust-lang.org/std/string/struct.String.html
    /// [`vcf::Header`]: header/struct.Header.html
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
        let mut header_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = match self.inner.fill_buf() {
                Ok(buf) => buf,
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            };

            if eol && !buf.is_empty() && buf[0] != HEADER_PREFIX {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == LINE_FEED) {
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

    /// Reads a single raw VCF record.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline character. The buffer can subsequently be parsed as a
    /// [`vcf::Record`].
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`records`]), but using this
    /// method allows control of the line buffer and whether the raw record should be parsed.
    ///
    /// If successful, the number of bytes is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// [`vcf::Record`]: record/struct.Record.html
    /// [`records`]: #method.records
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
    /// Unlike [`read_record`], each record is parsed as a [`vcf::Record`].
    ///
    /// [`read_record`]: #method.read_record
    /// [`vcf::Record`]: record/struct.Record.html
    ///
    /// # Examples
    ///
    /// ```no_run
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
    /// let mut records = reader.records();
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Reader<BufReader<bgzf::Reader<R>>>
where
    R: Read,
{
    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, BufReader};
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let data = Vec::new();
    /// let reader = vcf::Reader::new(BufReader::new(bgzf::Reader::new(&data[..])));
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(virtual_position.compressed(), 0);
    /// assert_eq!(virtual_position.uncompressed(), 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.get_ref().virtual_position()
    }
}

impl<R> Reader<BufReader<bgzf::Reader<R>>>
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
    /// # use std::{fs::File, io::{self, BufReader}};
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let mut reader = File::open("sample.vcf.gz")
    ///     .map(bgzf::Reader::new)
    ///     .map(BufReader::new)
    ///     .map(vcf::Reader::new)?;
    ///
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.get_mut().seek(pos)
    }
}

// Reads all bytes until a line feed ('\n') is reached.
//
// The buffer will not include the trailing newline ('\n' or '\r\n').
fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_line(buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            buf.pop();

            if buf.ends_with(CARRIAGE_RETURN) {
                buf.pop();
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
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
        let mut reader = Reader::new(DATA);

        let actual = reader.read_header()?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

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
        let mut reader = BufReader::new(&data[..]);
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        let data = b"noodles\r\n";
        let mut reader = BufReader::new(&data[..]);
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, "noodles");

        Ok(())
    }
}
