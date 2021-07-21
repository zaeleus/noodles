use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::{bin::Chunk, Metadata};

use super::{
    index::{reference_sequence, ReferenceSequence},
    Bin, Index, MAGIC_NUMBER,
};

/// A BAM index (BAI) reader.
///
/// A BAM index has three top-level fields:
///
///   1. a magic number,
///   2. a list of reference sequences,
///   3. and optionally, the number of unmapped reads in the associated BAM.
///
/// While these fields can be read individually, consider using [`super::read`] to read the entire
/// index at once.
///
/// # Examples
///
/// ```no_run
///# use std::{fs::File, io};
/// use noodles_bam::bai;
/// let mut reader = File::open("sample.bam.bai").map(bai::Reader::new)?;
/// reader.read_header()?;
/// let index = reader.read_index()?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Create a BAM index reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam::bai;
    /// let reader = File::open("sample.bam.bai").map(bai::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Read the BAM index header.
    ///
    /// The BAM index header is just the magic number of the file format.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam::bai;
    /// let mut reader = File::open("sample.bam.bai").map(bai::Reader::new)?;
    /// reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<()> {
        let mut magic = [0; 4];
        self.inner.read_exact(&mut magic)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAI header",
            ));
        }

        Ok(())
    }

    /// Reads the BAM index.
    ///
    /// The position of the stream is expected to be directly after the header.
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam::bai;
    /// let mut reader = File::open("sample.bam.bai").map(bai::Reader::new)?;
    /// reader.read_header()?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        let references = read_references(&mut self.inner)?;
        let n_no_coor = self.inner.read_u64::<LittleEndian>().ok();
        Ok(Index::new(references, n_no_coor))
    }
}

fn read_references<R>(reader: &mut R) -> io::Result<Vec<ReferenceSequence>>
where
    R: Read,
{
    let n_ref = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut references = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (bins, metadata) = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals, metadata));
    }

    Ok(references)
}

fn read_bins<R>(reader: &mut R) -> io::Result<(Vec<Bin>, Option<Metadata>)>
where
    R: Read,
{
    use reference_sequence::bin::METADATA_ID;

    let n_bin = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = Vec::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>()?;

        if id == METADATA_ID {
            metadata = read_metadata(reader).map(Some)?;
        } else {
            let chunks = read_chunks(reader)?;
            let bin = Bin::new(id, chunks);
            bins.push(bin);
        }
    }

    Ok((bins, metadata))
}

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

    for _ in 0..n_chunk {
        let chunk_beg = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        let chunk_end = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        chunks.push(Chunk::new(chunk_beg, chunk_end));
    }

    Ok(chunks)
}

fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: Read,
{
    let n_intv = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioffset = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioffset);
    }

    Ok(intervals)
}

fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: Read,
{
    use reference_sequence::bin::METADATA_CHUNK_COUNT;

    let n_chunk = reader.read_u32::<LittleEndian>()?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {}, got {}",
                METADATA_CHUNK_COUNT, n_chunk
            ),
        ));
    }

    let ref_beg = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let ref_end = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let n_mapped = reader.read_u64::<LittleEndian>()?;
    let n_unmapped = reader.read_u64::<LittleEndian>()?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() {
        let data = b"BAI\x01";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_header().is_ok());
    }

    #[test]
    fn test_read_header_with_invalid_magic_number() {
        let data = [];
        let mut reader = Reader::new(&data[..]);
        assert!(matches!(
            reader.read_header(),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAI";
        let mut reader = Reader::new(&data[..]);
        assert!(matches!(
            reader.read_header(),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = Reader::new(&data[..]);
        assert!(matches!(
            reader.read_header(),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_metadata() -> io::Result<()> {
        let data = [
            0x02, 0x00, 0x00, 0x00, // n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        let mut reader = &data[..];
        let actual = read_metadata(&mut reader)?;

        let expected = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
