use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_csi::{
    self as csi,
    binning_index::index::{
        reference_sequence::{index::LinearIndex, Bin, Metadata},
        ReferenceSequence,
    },
};

use super::{Index, MAGIC_NUMBER};

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
/// let index = reader.read_index()?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::bai;
    /// let reader = bai::Reader::new(io::empty());
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
    /// use noodles_bam::bai;
    /// let mut reader = bai::Reader::new(io::empty());
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
    /// use noodles_bam::bai;
    /// let reader = bai::Reader::new(io::empty());
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
    /// Creates a BAM index reader.
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

    /// Reads the BAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bam::bai;
    /// let mut reader = File::open("sample.bam.bai").map(bai::Reader::new)?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner)?;

        let references = read_references(&mut self.inner)?;
        let n_no_coor = read_unplaced_unmapped_record_count(&mut self.inner)?;

        let mut builder = Index::builder().set_reference_sequences(references);

        if let Some(n) = n_no_coor {
            builder = builder.set_unplaced_unmapped_record_count(n);
        }

        Ok(builder.build())
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
            "invalid BAI header",
        ))
    }
}

fn read_references<R>(reader: &mut R) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
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

fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    use csi::reader::index::reference_sequences::{bins::read_chunks, read_metadata};

    use super::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    fn duplicate_bin_error(id: usize) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)> {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("duplicate bin ID: {id}"),
        ))
    }

    let n_bin = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        if id == METADATA_ID {
            let m =
                read_metadata(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            if metadata.replace(m).is_some() {
                return duplicate_bin_error(id);
            }
        } else {
            let chunks =
                read_chunks(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let bin = Bin::new(chunks);

            if bins.insert(id, bin).is_some() {
                return duplicate_bin_error(id);
            }
        }
    }

    Ok((bins, metadata))
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

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: Read,
{
    match reader.read_u64::<LittleEndian>() {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic() {
        let data = b"BAI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());
    }

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAI";
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
    }

    #[test]
    fn test_read_bins() -> io::Result<()> {
        let data = [
            0x00, 0x00, 0x00, 0x00, // n_bin = 0
        ];
        let mut reader = &data[..];
        let (actual_bins, actual_metadata) = read_bins(&mut reader)?;
        assert!(actual_bins.is_empty());
        assert!(actual_metadata.is_none());

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        let mut reader = &data[..];
        let (actual_bins, actual_metadata) = read_bins(&mut reader)?;
        assert_eq!(actual_bins.len(), 1);
        assert!(actual_bins.get(&0).is_some());
        assert!(actual_metadata.is_some());

        let data = [
            0x01, 0x00, 0x00, 0x00, // n_bin = 1
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].n_chunk = 0
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x4a, 0x92, 0x00, 0x00, // bins[0].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[0].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_unmapped = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_unplaced_unmapped_record_count() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, None);

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, Some(8));

        Ok(())
    }
}
