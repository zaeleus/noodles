//! BAM reader.

mod builder;
pub mod header;
pub(crate) mod num;
pub(crate) mod query;
mod record;
mod record_buf;
mod record_bufs;
mod records;

use std::{
    ffi::CStr,
    io::{self, Read},
};

use bstr::BString;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{self as sam, alignment::RecordBuf, header::ReferenceSequences};

pub use self::{builder::Builder, query::Query, record_bufs::RecordBufs, records::Records};
use self::{record::read_record, record_buf::read_record_buf};
use crate::Record;

/// A BAM reader.
///
/// The BAM format is an encoded and compressed version of a SAM format.
///
/// The reader reads records sequentially but can use virtual positions to seek to offsets from the
/// start of a seekable stream.
///
/// # Examples
///
/// ## Read from a file
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_bam as bam;
///
/// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
/// let header = reader.read_header()?;
///
/// for result in reader.records() {
///     let record = result?;
///     // ...
/// }
/// # Ok::<_, io::Error>(())
/// ```
///
/// ## Use a custom BGZF decoder
///
/// [`Reader::new`] wraps the input stream with a default BGZF decoder. This can be swapped for a
/// custom decoder, e.g., [`flate2::read::MultiGzDecoder`],
/// [`noodles_bgzf::io::MultithreadedReader`], etc., using [`Reader::from`].
///
/// [`flate2::read::MultiGzDecoder`]: https://docs.rs/flate2/latest/flate2/read/struct.MultiGzDecoder.html
///
/// ### `flate2::read::MultiGzDecoder`
///
/// ```
/// # use std::{fs::File, io};
/// use flate2::read::MultiGzDecoder;
/// use noodles_bam as bam;
///
/// let decoder = MultiGzDecoder::new(io::empty());
/// let _reader = bam::io::Reader::from(decoder);
/// ```
///
/// ### `noodles_bgzf::io::MultithreadedReader`
///
/// ```
/// # use std::{fs::File, io, num::NonZero, thread};
/// use noodles_bam as bam;
/// use noodles_bgzf as bgzf;
///
/// let worker_count = thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
/// let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, io::empty());
/// let _reader = bam::io::Reader::from(decoder);
/// ```
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
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let reader = bam::io::Reader::from(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let mut reader = bam::io::Reader::from(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let reader = bam::io::Reader::from(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Returns a BAM header reader.
    ///
    /// This creates an adapter that reads at most the length of the header, i.e., the BAM magic
    /// number, the SAM header, and reference sequences.
    ///
    /// It is more ergonomic to read the BAM header as a SAM header using [`Self::read_header`],
    /// but this adapter allows for control of how the header is read, e.g., to read the raw SAM
    /// header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::Read};
    /// use noodles_bam as bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    ///
    /// let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;
    /// let mut raw_header = String::new();
    /// raw_sam_header_reader.read_to_string(&mut raw_header)?;
    /// raw_sam_header_reader.discard_to_end()?;
    ///
    /// header_reader.read_reference_sequences()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
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

    /// Reads a record into an alignment record buffer.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into an internal buffer. This buffer is then used to decode fields into the given record.
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
        _header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, &mut self.buf, record)
    }

    /// Reads a record.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into the record's buffer. No fields are decoded, meaning the record is not necessarily
    /// valid. However, the structure of the buffer is guaranteed to be record-like.
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
        let fields = record.fields_mut();

        let block_size = match read_record(&mut self.inner, &mut fields.buf)? {
            0 => return Ok(0),
            n => n,
        };

        fields.index()?;

        Ok(block_size)
    }

    /// Returns an iterator over alignment record buffers starting from the current stream
    /// position.
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
    pub fn record_bufs<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> RecordBufs<'r, 'h, R> {
        RecordBufs::new(self, header)
    }

    /// Returns an iterator over records.
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

impl<R> Reader<bgzf::io::Reader<R>>
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
        Self::from(bgzf::io::Reader::new(reader))
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

    /// Returns a reader over records that intersect the given region.
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
    /// let index = bai::fs::read("sample.bam.bai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, I>(
        &'r mut self,
        header: &sam::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'r, R>>
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
    /// let index = bai::fs::read("sample.bam.bai")?;
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
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + use<'r, I, R>>
    where
        I: BinningIndex,
    {
        if let Some(pos) = index.last_first_record_start_position() {
            self.get_mut().seek_to_virtual_position(pos)?;
        } else {
            self.seek_to_first_record()?;
        }

        Ok(self.records().filter(|result| {
            result
                .as_ref()
                .map(|record| record.flags().is_unmapped())
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
        _header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a> {
        Box::new(
            self.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            }),
        )
    }
}

pub(crate) fn bytes_with_nul_to_bstring(buf: &[u8]) -> io::Result<BString> {
    CStr::from_bytes_with_nul(buf)
        .map(|c_str| c_str.to_bytes().into())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
