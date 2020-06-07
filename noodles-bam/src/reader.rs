//! BAM reading and iterators

mod query;
mod records;

pub use self::{query::Query, records::Records};

use std::{
    ffi::CStr,
    io::{self, Read, Seek},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles::Region;
use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_sam::header::{ReferenceSequence, ReferenceSequences};

use super::{bai, Record, MAGIC_NUMBER};

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
/// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
/// reader.read_header()?;
/// reader.read_reference_sequences()?;
///
/// for result in reader.records() {
///     let record = result?;
///     println!("{:?}", record);
/// }
///
/// # Ok::<(), io::Error>(())
/// ```
pub struct Reader<R>
where
    R: Read,
{
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BAM reader.
    ///
    /// The given reader must be a raw BGZF stream, as the underlying reader wraps it in a decoder.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// let mut _reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    /// Reads the raw SAM header.
    ///
    /// The BAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a `String`. It can subsequently be parsed as a
    /// `sam::Header`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        let magic = read_magic(&mut self.inner)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM header",
            ));
        }

        read_header(&mut self.inner)
    }

    /// Reads the binary reference sequences after the SAM header.
    ///
    /// This is not the same as the `@SQ` records in the SAM header. A BAM has a list of reference
    /// sequences containing name and length tuples after the SAM header and before the list of
    /// records.
    ///
    /// The position of the stream is expected to be directly after the header.
    ///
    /// This returns a list of `sam::header::ReferenceSequence` objects, which can be used to build
    /// a minimal `sam::Header` if the SAM header is empty.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// let reference_sequences = reader.read_reference_sequences()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_reference_sequences(&mut self) -> io::Result<Vec<ReferenceSequence>> {
        let n_ref = self.inner.read_u32::<LittleEndian>()?;
        let mut reference_sequences = Vec::with_capacity(n_ref as usize);

        for _ in 0..n_ref {
            let reference_sequence = read_reference_sequence(&mut self.inner)?;
            reference_sequences.push(reference_sequence);
        }

        Ok(reference_sequences)
    }

    /// Reads a single record.
    ///
    /// The record block size (`bs`) is read from the underlying stream, and `bs` addition bytes
    /// are read into the given record.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// It is more ergonomic to read records using an iterator (see `records` and `query`), but
    /// using this method directly allows the reuse of a single `bam::Record` buffer.
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
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// reader.read_reference_sequences()?;
    ///
    /// let mut record = bam::Record::default();
    /// reader.read_record(&mut record)?;
    ///
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, mut record: &mut Record) -> io::Result<usize> {
        let block_size = match self.inner.read_u32::<LittleEndian>() {
            Ok(bs) => bs as usize,
            Err(e) => match e.kind() {
                io::ErrorKind::UnexpectedEof => return Ok(0),
                _ => return Err(e),
            },
        };

        record.resize(block_size);
        self.inner.read_exact(&mut record)?;

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
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// reader.read_reference_sequences()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
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
    /// let reader = bam::Reader::new(&data[..]);
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(virtual_position.compressed(), 0);
    /// assert_eq!(virtual_position.uncompressed(), 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the underlying BGZF reader to the given virtual position.
    ///
    /// Virtual positions typically come from the associated BAM index file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    ///
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        self.inner.seek(pos)
    }

    /// Returns an iterator over records that intersect the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles::Region;
    /// use noodles_bam::{self as bam, bai};
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// let header: sam::Header = reader.read_header()?.parse()?;
    ///
    /// let reference_sequences = header.reference_sequences();
    /// let index = bai::read("sample.bam.bai")?;
    /// let region = Region::mapped("sq0", 17711, Some(28657));
    /// let query = reader.query(&reference_sequences, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(
        &mut self,
        reference_sequences: &ReferenceSequences,
        index: &bai::Index,
        region: &Region,
    ) -> io::Result<Query<R>> {
        let (i, _, start, end) = region.resolve(reference_sequences).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("could not resolve region: {:?}", region),
            )
        })?;

        let index_reference_sequence = index.reference_sequences().get(i).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "could not find reference in index: {} >= {}",
                    i,
                    reference_sequences.len()
                ),
            )
        })?;

        let query_bins = index_reference_sequence.query(start, end);

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .cloned()
            .collect();

        let min_offset = index_reference_sequence.min_offset(start);
        let merged_chunks = bai::optimize_chunks(&chunks, min_offset);

        Ok(Query::new(self, merged_chunks, i, start, end))
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    Ok(magic)
}

fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_u32::<LittleEndian>()?;

    let mut c_text = vec![0; l_text as usize];
    reader.read_exact(&mut c_text)?;

    // Headers are not necessarily NUL-terminated.
    bytes_with_nul_to_string(&c_text).or_else(|_| {
        String::from_utf8(c_text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence>
where
    R: Read,
{
    let l_name = reader.read_u32::<LittleEndian>()?;

    let mut c_name = vec![0; l_name as usize];
    reader.read_exact(&mut c_name)?;

    let name = bytes_with_nul_to_string(&c_name)?;
    let l_ref = reader.read_u32::<LittleEndian>()?;

    Ok(ReferenceSequence::new(name, l_ref as i32))
}

fn bytes_with_nul_to_string(buf: &[u8]) -> io::Result<String> {
    CStr::from_bytes_with_nul(buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_str| {
            c_str
                .to_str()
                .map(|s| s.to_string())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}
