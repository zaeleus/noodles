//! BAM reader and iterators.

mod lazy_records;
pub(crate) mod query;
pub mod record;
mod records;
mod unmapped_records;

pub use self::{
    lazy_records::LazyRecords, query::Query, records::Records, unmapped_records::UnmappedRecords,
};

use std::{
    ffi::CStr,
    io::{self, Read, Seek},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::Record,
    header::{
        record::value::{map::ReferenceSequence, Map},
        ReferenceSequences,
    },
};

use super::{bai, lazy, MAGIC_NUMBER};

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
    /// let reader = bam::Reader::from(&data[..]);
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
    /// let mut reader = bam::Reader::from(&data[..]);
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
    /// let reader = bam::Reader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the raw SAM header.
    ///
    /// The BAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`noodles_sam::Header`].
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
        read_magic(&mut self.inner)?;
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
    /// This returns a reference sequence dictionary ([`noodles_sam::header::ReferenceSequences`]),
    /// which can be used to build a minimal [`noodles_sam::Header`] if the SAM header is empty.
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
    pub fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(&mut self.inner)
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
    /// [`Self::query`]), but using this method directly allows the reuse of a single [`Record`]
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
    /// use noodles_sam::alignment::Record;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// reader.read_reference_sequences()?;
    ///
    /// let mut record = Record::default();
    /// reader.read_record(&mut record)?;
    ///
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        use self::record::read_record;
        read_record(&mut self.inner, &mut self.buf, record)
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
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// reader.read_reference_sequences()?;
    ///
    /// let mut record = bam::lazy::Record::default();
    /// reader.read_lazy_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        use self::record::read_block_size;

        let block_size = match read_block_size(&mut self.inner)? {
            0 => return Ok(0),
            n => n,
        };

        record.buf.resize(block_size, 0);
        self.inner.read_exact(&mut record.buf)?;

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
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
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
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// reader.read_header()?;
    /// reader.read_reference_sequences()?;
    ///
    /// for result in reader.lazy_records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn lazy_records(&mut self) -> LazyRecords<'_, R> {
        LazyRecords::new(self)
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
    /// let reader = bam::Reader::new(&data[..]);
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
    /// let reader = bam::Reader::new(&data[..]);
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
    /// let mut reader = bam::Reader::new(Cursor::new(Vec::new()));
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }

    // Seeks to the first record by setting the cursor to the beginning of the stream and
    // (re)reading the header and binary reference sequences.
    fn seek_to_first_record(&mut self) -> io::Result<bgzf::VirtualPosition> {
        self.seek(bgzf::VirtualPosition::default())?;
        self.read_header()?;
        self.read_reference_sequences()?;
        Ok(self.virtual_position())
    }

    /// Returns an iterator over records that intersect the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bam::{self as bam, bai};
    /// use noodles_core::Region;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// let header: sam::Header = reader.read_header()?.parse()?;
    ///
    /// let reference_sequences = header.reference_sequences();
    /// let index = bai::read("sample.bam.bai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(reference_sequences, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<I>(
        &mut self,
        reference_sequences: &ReferenceSequences,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'_, R>>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(reference_sequences, region)?;

        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self,
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
    /// let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
    /// let index = bai::read("sample.bam.bai")?;
    /// let query = reader.query_unmapped(&index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn query_unmapped(&mut self, index: &bai::Index) -> io::Result<UnmappedRecords<'_, R>> {
        if let Some(pos) = index.first_record_in_last_linear_bin_start_position() {
            self.seek(pos)?;
        } else {
            self.seek_to_first_record()?;
        }

        Ok(UnmappedRecords::new(self))
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

impl<R> sam::AlignmentReader<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        read_alignment_header(&mut self.inner)
    }

    fn alignment_records<'a>(
        &'a mut self,
        _: &'a fasta::Repository,
        _: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Record>> + 'a> {
        Box::new(self.records())
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut text = vec![0; l_text];
    reader.read_exact(&mut text)?;

    // § 4.2 The BAM format (2021-06-03): "Plain header text in SAM; not necessarily
    // NUL-terminated".
    bytes_with_nul_to_string(&text).or_else(|_| {
        String::from_utf8(text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: Read,
{
    let n_ref = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader)?;
        let name = reference_sequence.name().to_string();
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_reference_sequence<R>(reader: &mut R) -> io::Result<Map<ReferenceSequence>>
where
    R: Read,
{
    let l_name = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name)?;

    let name = bytes_with_nul_to_string(&c_name).and_then(|name| {
        name.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let l_ref = reader.read_u32::<LittleEndian>().and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Map::<ReferenceSequence>::new(name, l_ref)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_alignment_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    read_magic(reader)?;

    let header = read_header(reader).and_then(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    read_reference_sequences(reader)?;

    Ok(header)
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
                    "region reference sequence does not exist in reference sequences: {:?}",
                    region
                ),
            )
        })
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

    use super::*;

    #[test]
    fn test_read_magic() -> io::Result<()> {
        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
        let expected = "@HD\tVN:1.6\n";

        let data_len = expected.len() as u32;
        let mut data = data_len.to_le_bytes().to_vec();
        data.extend(expected.as_bytes());

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::record::value::map::reference_sequence::Name;

        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader)?;

        let expected: ReferenceSequences = [("sq0".parse()?, 8)]
            .into_iter()
            .map(|(name, len): (Name, usize)| {
                let sn = name.to_string();
                Map::<ReferenceSequence>::new(name, len).map(|rs| (sn, rs))
            })
            .collect::<Result<_, _>>()?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_alignment_header() -> Result<(), Box<dyn std::error::Error>> {
        use bytes::BufMut;
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        let mut data = Vec::new();
        data.put_slice(MAGIC_NUMBER); // magic
        data.put_u32_le(11); // l_text
        data.put_slice(b"@HD\tVN:1.6\n"); // text
        data.put_u32_le(1); // n_ref
        data.put_u32_le(4); // ref[0].l_name
        data.put_slice(b"sq0\x00"); // ref[0].name
        data.put_u32_le(8); // ref[0].l_ref

        let mut reader = &data[..];
        let actual = read_alignment_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
