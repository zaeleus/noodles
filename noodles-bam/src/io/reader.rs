//! BAM reader.

mod builder;
mod header;
mod query;
mod record;
mod record_buf;
mod record_bufs;
mod records;

use std::{
    ffi::CStr,
    io::{self, Read, Seek},
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{self as sam, alignment::RecordBuf, header::ReferenceSequences};

pub use self::{builder::Builder, query::Query, record_bufs::RecordBufs, records::Records};
use self::{record::read_record, record_buf::read_record_buf};
use crate::Record;

/// A BAM reader.
///
/// A BAM file is an encoded and compressed version of a SAM file. While a SAM file has a header
/// and a list of records, a BAM is comprised of three parts:
///
///   1. a SAM header,
///   2. a list of reference sequences, and
///   3. a list of encoded SAM records.
///
/// The reader reads records sequentially but can use virtual positions to seek to offsets from the
/// start of a seekable stream.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_bam as bam;
///
/// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
/// let header = reader.read_header()?;
///
/// for result in reader.record_bufs(&header) {
///     let record = result?;
///     // ...
/// }
/// # Ok::<_, io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::io::Reader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let mut reader = bam::io::Reader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::io::Reader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the SAM header.
    ///
    /// This verifies the BAM magic number, reads and parses the raw SAM header, and reads the
    /// binary reference sequences. If the SAM header has a reference sequence dictionary, it must
    /// match the binary reference sequences; otherwise, the binary reference sequences are added
    /// to the SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        use self::header::read_header;
        read_header(&mut self.inner)
    }

    /// Reads a single record.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into an internal buffer. This buffer is used to populate the given record.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`] and
    /// [`Self::query`]), but using this method directly allows the reuse of a single [`RecordBuf`]
    /// buffer.
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// use noodles_sam::alignment::RecordBuf;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let mut record = RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record)?;
    ///
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, header, &mut self.buf, record)
    }

    /// Reads a single record without eagerly decoding its fields.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into the lazy record's buffer. No fields are decoded, meaning the record is not necessarily
    /// valid. However, the structure of the byte stream is guaranteed to be record-like.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut record = bam::Record::default();
    /// reader.read_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let block_size = match read_record(&mut self.inner, &mut record.buf)? {
            0 => return Ok(0),
            n => n,
        };

        record.index()?;

        Ok(block_size)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.record_bufs(&header) {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn record_bufs<'a>(&'a mut self, header: &'a sam::Header) -> RecordBufs<'_, R> {
        RecordBufs::new(self, header)
    }

    /// Returns an iterator over lazy records.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Creates a BAM reader.
    ///
    /// The given reader must be a raw BGZF stream, as the underlying reader wraps it in a decoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::io::Reader::new(&data[..]);
    /// ```
    pub fn new(reader: R) -> Self {
        Self::from(bgzf::Reader::new(reader))
    }

    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    ///
    /// let data = Vec::new();
    /// let reader = bam::io::Reader::new(&data[..]);
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
    /// Seeks the underlying BGZF reader to the given virtual position.
    ///
    /// Virtual positions typically come from the associated BAM index file.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    /// let mut reader = bam::io::Reader::new(Cursor::new(Vec::new()));
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }

    // Seeks to the first record by setting the cursor to the beginning of the stream and
    // (re)reading the header.
    fn seek_to_first_record(&mut self) -> io::Result<bgzf::VirtualPosition> {
        self.seek(bgzf::VirtualPosition::default())?;
        self.read_header()?;
        Ok(self.virtual_position())
    }

    /// Returns an iterator over records that intersect the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bam::{self as bam, bai};
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let index = bai::read("sample.bam.bai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<I>(
        &mut self,
        header: &sam::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'_, R>>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(header.reference_sequences(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self.get_mut(),
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam::{self as bam, bai};
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let index = bai::read("sample.bam.bai")?;
    /// let query = reader.query_unmapped(&index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn query_unmapped<'r, I>(
        &'r mut self,
        index: &I,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'r>
    where
        I: BinningIndex,
    {
        use sam::alignment::record::Flags;

        if let Some(pos) = index.last_first_record_start_position() {
            self.seek(pos)?;
        } else {
            self.seek_to_first_record()?;
        }

        Ok(self.records().filter(|result| {
            result
                .as_ref()
                .map(|record| Flags::from(record.flags()).is_unmapped())
                .unwrap_or(true)
        }))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}

impl<R> sam::alignment::io::Read<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        self.read_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a> {
        Box::new(
            self.record_bufs(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            }),
        )
    }
}

pub(crate) fn bytes_with_nul_to_string(buf: &[u8]) -> io::Result<String> {
    CStr::from_bytes_with_nul(buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_str| {
            c_str
                .to_str()
                .map(|s| s.to_string())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
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
