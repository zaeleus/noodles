//! SAM reader.

mod builder;
pub mod header;
pub(crate) mod query;
mod record;
pub(crate) mod record_buf;
mod record_bufs;

use std::{
    io::{self, BufRead},
    iter,
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;

pub(crate) use self::record::read_record;
pub use self::{builder::Builder, record_bufs::RecordBufs};
use self::{header::read_header, query::Query, record_buf::read_record_buf};
use crate::{Header, Record, alignment::RecordBuf, header::ReferenceSequences};

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
/// use noodles_sam as sam;
///
/// let mut reader = sam::io::reader::Builder::default().build_from_path("sample.sam")?;
/// let header = reader.read_header()?;
///
/// for result in reader.records() {
///     let record = result?;
///     // ...
/// }
/// # Ok::<_, std::io::Error>(())
/// ```
#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::io::Reader::new(&data[..]);
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
    /// let mut reader = sam::io::Reader::new(&data[..]);
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
    /// let reader = sam::io::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
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
    /// # use std::io;
    /// use noodles_sam as sam;
    /// let reader = sam::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self::from(inner)
    }

    /// Returns a SAM header reader.
    ///
    /// This creates an adapter that reads at most the length of the header, i.e., all lines
    /// prefixed with a `@` (at sign).
    ///
    /// It is more ergonomic to read and parse the header using [`Self::read_header`], but using
    /// this adapter allows for control of how the header is read, e.g., to read the raw SAM
    /// header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::Read;
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// let mut header_reader = reader.header_reader();
    ///
    /// let mut raw_header = String::new();
    /// header_reader.read_to_string(&mut raw_header)?;
    ///
    /// assert_eq!(raw_header, "@HD\tVN:1.6\n");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// The SAM header is optional, and if it is missing, an empty [`Header`] is returned.
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
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// let actual = reader.read_header()?;
    ///
    /// let expected = sam::Header::builder()
    ///     .set_header(Default::default())
    ///     .build();
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<Header> {
        read_header(&mut self.inner)
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// This reads a line from the underlying stream until a newline is reached and parses that
    /// line into the given record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`] and
    /// [`Self::query`]), but using this method directly allows reuse of a [`RecordBuf`].
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut record = RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record)?;
    ///
    /// assert_eq!(record, RecordBuf::default());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_record_buf(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, &mut self.buf, header, record)
    }

    /// Returns an iterator over alignment record buffers starting from the current stream
    /// position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
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
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut records = reader.record_bufs(&header);
    /// assert!(records.next().is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn record_bufs<'a>(&'a mut self, header: &'a Header) -> RecordBufs<'a, R> {
        RecordBufs::new(self, header)
    }

    /// Reads a record.
    ///
    /// This reads SAM fields from the underlying stream into the given record's buffer until a
    /// newline is reached. No fields are decoded, meaning the record is not necessarily valid.
    /// However, the structure of the buffer is guaranteed to be record-like.
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
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut record = sam::Record::default();
    /// reader.read_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record)
    }

    /// Returns an iterator over records.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
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
    /// let mut reader = sam::io::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn records(&mut self) -> impl Iterator<Item = io::Result<Record>> {
        let mut record = Record::default();

        iter::from_fn(move || match self.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record.clone())),
            Err(e) => Some(Err(e)),
        })
    }
}

impl<R> Reader<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    // Seeks to the first record by setting the cursor to the beginning of the stream and
    // (re)reading the header.
    fn seek_to_first_record(&mut self) -> io::Result<bgzf::VirtualPosition> {
        self.get_mut()
            .seek_to_virtual_position(bgzf::VirtualPosition::default())?;

        self.read_header()?;

        Ok(self.get_ref().virtual_position())
    }

    /// Returns an iterator over records that intersect the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.sam.gz")
    ///     .map(bgzf::io::Reader::new)
    ///     .map(sam::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    ///
    /// let index = csi::fs::read("sample.sam.gz.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r, I>(
        &'r mut self,
        header: &'h Header,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + use<'r, 'h, I, R>>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(header.reference_sequences(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self.get_mut(),
            chunks,
            header,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.sam.gz")
    ///     .map(bgzf::io::Reader::new)
    ///     .map(sam::io::Reader::new)?;
    ///
    /// reader.read_header()?;
    ///
    /// let index = csi::fs::read("sample.sam.gz.csi")?;
    /// let query = reader.query_unmapped(&index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn query_unmapped<'r, I>(
        &'r mut self,
        index: &I,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + use<'r, I, R>>
    where
        I: BinningIndex,
    {
        if let Some(pos) = index.last_first_record_start_position() {
            self.get_mut().seek_to_virtual_position(pos)?;
        } else {
            self.seek_to_first_record()?;
        }

        let mut record = Record::default();

        Ok(iter::from_fn(move || {
            loop {
                match self.read_record(&mut record) {
                    Ok(0) => return None,
                    Ok(_) => {
                        let result = record.flags().map(|flags| flags.is_unmapped());

                        match result {
                            Ok(true) => return Some(Ok(record.clone())),
                            Ok(false) => {}
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
        }))
    }
}

impl<R> From<R> for Reader<R>
where
    R: BufRead,
{
    fn from(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}

impl<R> crate::alignment::io::Read<R> for Reader<R>
where
    R: BufRead,
{
    fn read_alignment_header(&mut self) -> io::Result<Header> {
        self.read_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        _header: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn crate::alignment::Record>>> + 'a> {
        Box::new(self.records().map(|result| {
            result.map(|record| Box::new(record) as Box<dyn crate::alignment::Record>)
        }))
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
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

pub(crate) fn resolve_region(
    reference_sequences: &ReferenceSequences,
    region: &Region,
) -> io::Result<usize> {
    reference_sequences
        .get_index_of(region.name())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "region reference sequence does not exist in reference sequences: {region:?}"
                ),
            )
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_line() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut reader: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_line(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"noodles\n", b"noodles")?;
        t(&mut buf, b"noodles\r\n", b"noodles")?;
        t(&mut buf, b"noodles", b"noodles")?;

        Ok(())
    }
}
