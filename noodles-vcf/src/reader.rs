//! VCF reader and iterators.

mod builder;
mod header;
pub(crate) mod query;
pub mod record;
mod records;

use crate::lazy;

pub(crate) use self::record::parse_record;
pub use self::{builder::Builder, query::Query, records::Records};

use std::{
    io::{self, BufRead, Read, Seek},
    str,
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;

use self::header::read_header;
use super::{Header, Record, VariantReader};

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
/// use noodles_vcf as vcf;
///
/// let mut reader = vcf::reader::Builder::default().build_from_path("sample.vcf")?;
/// let header = reader.read_header()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     println!("{:?}", record);
/// }
/// # Ok::<_, std::io::Error>(())
/// ```
#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
    buf: String,
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
        Self {
            inner,
            buf: String::new(),
        }
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

    /// Reads the VCF header.
    ///
    /// This reads all header lines prefixed with a `#` (number sign), which includes the header
    /// header (`#CHROM`...), and parses it as a [`crate::Header`].
    ///
    /// The position of the stream is expected to be at the start.
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
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<Header> {
        read_header(&mut self.inner)
    }

    /// Reads a single VCF record.
    ///
    /// This reads a line from the underlying stream until a newline is reached and parses that
    /// line into the given record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`]), but using
    /// this method allows control of the record buffer.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
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
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut record = vcf::Record::default();
    /// reader.read_record(&header, &mut record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_record(&mut self, header: &Header, record: &mut Record) -> io::Result<usize> {
        self.buf.clear();

        match read_line(&mut self.inner, &mut self.buf)? {
            0 => Ok(0),
            n => {
                parse_record(&self.buf, header, record)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                Ok(n)
            }
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
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
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut records = reader.records(&header);
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(&'r mut self, header: &'h Header) -> Records<'r, 'h, R> {
        Records::new(self, header)
    }

    /// Reads a single record without eagerly parsing its fields.
    ///
    /// The reads VCF record fields from the underlying stream into the given record's buffer until
    /// a newline is reached. No fields are parsed, meaning the record is not necessarily valid.
    /// However, the structure of the line is guaranteed to be record-like.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
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
    /// let mut reader = vcf::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut record = vcf::lazy::Record::default();
    /// reader.read_lazy_record(&mut record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        read_lazy_record(&mut self.inner, record)
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
    /// let header = reader.read_header()?;
    ///
    /// let index = tabix::read("sample.vcf.gz.tbi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h, I>(
        &'r mut self,
        header: &'h Header,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, R>>
    where
        I: BinningIndex,
    {
        let (reference_sequence_id, reference_sequence_name) = resolve_region(index, region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self.get_mut(),
            chunks,
            reference_sequence_name,
            region.interval(),
            header,
        ))
    }
}

impl<R> VariantReader<R> for Reader<R>
where
    R: BufRead,
{
    fn read_variant_header(&mut self) -> io::Result<Header> {
        self.read_header()
    }

    fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<crate::Record>> + 'r> {
        Box::new(self.records(header))
    }
}

// Reads all bytes until a line feed ('\n') or EOF is reached.
//
// The buffer will not include the trailing newline ('\n' or '\r\n').
fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

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

fn read_lazy_record<R>(reader: &mut R, record: &mut lazy::Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.buf.clear();

    let mut len = 0;

    len += read_field(reader, &mut record.buf)?;
    record.bounds.chromosome_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.position_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.ids_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.reference_bases_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.alternate_bases_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.quality_score_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.filters_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.info_end = record.buf.len();

    len += read_line(reader, &mut record.buf)?;

    Ok(len)
}

fn read_field<R>(reader: &mut R, dst: &mut String) -> io::Result<usize>
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

        let (buf, n) = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                is_delimiter = true;
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst.push_str(s);

        len += n;

        reader.consume(n);
    }

    Ok(len)
}

pub(crate) fn resolve_region<I>(index: &I, region: &Region) -> io::Result<(usize, Vec<u8>)>
where
    I: BinningIndex,
{
    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tabix header"))?;

    let region_name = str::from_utf8(region.name())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let i = header
        .reference_sequence_names()
        .get_index_of(region_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "region reference sequence does not exist in reference sequences: {region:?}"
                ),
            )
        })?;

    Ok((i, region.name().into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        static DATA: &[u8] = b"\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t1\t.\tA\t.\t.\tPASS\t.
";

        let mut reader = Reader::new(DATA);
        let header = reader.read_header()?;

        let mut record = Record::default();

        let bytes_read = reader.read_record(&header, &mut record)?;
        assert_eq!(bytes_read, 21);

        let bytes_read = reader.read_record(&header, &mut record)?;
        assert_eq!(bytes_read, 0);

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

    #[test]
    fn test_read_lazy_record() -> io::Result<()> {
        let mut src = &b"sq0\t1\t.\tA\t.\t.\tPASS\t."[..];

        let mut record = lazy::Record::default();
        read_lazy_record(&mut src, &mut record)?;

        assert_eq!(record.buf, "sq01.A..PASS.");

        assert_eq!(record.bounds.chromosome_end, 3);
        assert_eq!(record.bounds.position_end, 4);
        assert_eq!(record.bounds.ids_end, 5);
        assert_eq!(record.bounds.reference_bases_end, 6);
        assert_eq!(record.bounds.alternate_bases_end, 7);
        assert_eq!(record.bounds.quality_score_end, 8);
        assert_eq!(record.bounds.filters_end, 12);
        assert_eq!(record.bounds.info_end, 13);

        Ok(())
    }
}
