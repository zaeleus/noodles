use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use super::{
    index::{
        reference_sequence::{self, Metadata},
        ReferenceSequence,
    },
    Bin, Chunk, Index, MAGIC_NUMBER,
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
    let n_ref = reader.read_u32::<LittleEndian>()?;
    let mut references = Vec::with_capacity(n_ref as usize);

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
    use reference_sequence::metadata::TryFromBinError;

    let n_bin = reader.read_u32::<LittleEndian>()?;

    let mut bins = Vec::with_capacity(n_bin as usize);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>()?;
        let chunks = read_chunks(reader)?;

        let bin = Bin::new(id, chunks);

        match Metadata::try_from(&bin) {
            Ok(m) => metadata = Some(m),
            Err(TryFromBinError::InvalidMagicNumber(_)) => bins.push(bin),
            Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
        }
    }

    Ok((bins, metadata))
}

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_u32::<LittleEndian>()?;
    let mut chunks = Vec::with_capacity(n_chunk as usize);

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
    let n_intv = reader.read_u32::<LittleEndian>()?;
    let mut intervals = Vec::with_capacity(n_intv as usize);

    for _ in 0..n_intv {
        let ioffset = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioffset);
    }

    Ok(intervals)
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
        assert!(reader.read_header().is_err());

        let data = b"BAI";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_header().is_err());

        let data = b"MThd";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_header().is_err());
    }
}
