//! BCF reader.

mod builder;
pub mod header;
pub(crate) mod num;
pub(crate) mod query;
pub(crate) mod record;
pub(crate) mod record_buf;
mod record_bufs;

pub use self::{builder::Builder, query::Query, record_bufs::RecordBufs};

use std::{
    io::{self, BufRead, Read},
    iter, str,
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_vcf::{self as vcf, header::string_maps::ContigStringMap, variant::RecordBuf};

use self::{header::read_header, record::read_record, record_buf::read_record_buf};
use crate::Record;

/// A BCF reader.
///
/// The BCF format is comprised of two parts: 1) a VCF header and 2) a list of records.
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
    /// use noodles_bcf as bcf;
    /// let reader = bcf::io::Reader::from(io::empty());
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
    /// use noodles_bcf as bcf;
    /// let mut reader = bcf::io::Reader::from(io::empty());
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
    /// use noodles_bcf as bcf;
    /// let reader = bcf::io::Reader::from(io::empty());
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
    /// Returns a BCF header reader.
    ///
    /// This creates an adapter that reads at most the length of the header, i.e., the BCF magic
    /// number, the format version, and VCF header.
    ///
    /// It is more ergonomic to read the BCF header as a VCF header using [`Self::read_header`],
    /// but this adapter allows for control of how the header is read, e.g., to read the raw VCF
    /// header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::Read};
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    ///
    /// let mut raw_vcf_header_reader = header_reader.raw_vcf_header_reader()?;
    /// let mut raw_header = String::new();
    /// raw_vcf_header_reader.read_to_string(&mut raw_header)?;
    /// raw_vcf_header_reader.discard_to_end()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the VCF header.
    ///
    /// This verifies the BCF magic number, discards the file format version, and reads and parses
    /// the raw VCF header. Associated string maps are also built from the raw header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        read_header(&mut self.inner)
    }

    /// Reads a single record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`]), but using
    /// this method directly allows the reuse of a single [`vcf::Record`] buffer.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record_buf(
        &mut self,
        header: &vcf::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, header, &mut self.buf, record)
    }

    /// Reads a single record without eagerly decoding (most of) its fields.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergnomic to read records using an iterator (see [`Self::records`]), but using
    /// this method directly allows the reuse of a single [`Record`] buffer.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut record = bcf::Record::default();
    /// reader.read_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record)
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
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.record_bufs(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    pub fn record_bufs<'r, 'h>(&'r mut self, header: &'h vcf::Header) -> RecordBufs<'r, 'h, R> {
        RecordBufs::new(self, header)
    }

    /// Returns an iterator over lazy records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
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

impl<R> Reader<bgzf::io::Reader<R>>
where
    R: Read,
{
    /// Creates a BCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let reader = bcf::io::Reader::new(io::empty());
    /// ```
    pub fn new(reader: R) -> Self {
        Self::from(bgzf::io::Reader::new(reader))
    }
}

impl<R> Reader<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Returns a reader over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    /// use noodles_core::Region;
    /// use noodles_csi as csi;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let index = csi::fs::read("sample.bcf.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h, I>(
        &'r mut self,
        header: &'h vcf::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, R>>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(header.string_maps().contigs(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            &mut self.inner,
            header,
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
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

impl<R> vcf::variant::io::Read<R> for Reader<R>
where
    R: BufRead,
{
    fn read_variant_header(&mut self) -> io::Result<vcf::Header> {
        self.read_header()
    }

    fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        _: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn vcf::variant::Record>>> + 'r> {
        Box::new(
            self.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            }),
        )
    }
}

pub(crate) fn resolve_region(
    contig_string_map: &ContigStringMap,
    region: &Region,
) -> io::Result<usize> {
    let region_name = str::from_utf8(region.name())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    contig_string_map.get_index_of(region_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("region does not exist in contigs: {region:?}"),
        )
    })
}
