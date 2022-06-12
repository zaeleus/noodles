use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;

use super::{
    binning_index::ReferenceSequenceExt,
    index::{
        reference_sequence::{bin::Chunk, Bin, Metadata},
        ReferenceSequence,
    },
    BinningIndex, Index, MAGIC_NUMBER,
};

/// A CSI writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CSI writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let writer = csi::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Writes a coordinate-sorted index (CSI).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let mut writer = csi::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_magic(&mut self.inner)?;

        let min_shift = i32::from(index.min_shift());
        self.inner.write_i32::<LittleEndian>(min_shift)?;

        let depth = i32::from(index.depth());
        self.inner.write_i32::<LittleEndian>(depth)?;

        write_aux(&mut self.inner, index.aux())?;
        write_reference_sequences(&mut self.inner, index.depth(), index.reference_sequences())?;

        if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
            self.inner.write_u64::<LittleEndian>(n_no_coor)?;
        }

        Ok(())
    }
}

fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)
}

fn write_aux<W>(writer: &mut W, aux: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let l_aux =
        i32::try_from(aux.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(l_aux)?;

    writer.write_all(aux)?;

    Ok(())
}

fn write_reference_sequences<W>(
    writer: &mut W,
    depth: u8,
    reference_sequences: &[ReferenceSequence],
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_ref)?;

    for reference_sequence in reference_sequences {
        write_bins(
            writer,
            depth,
            reference_sequence.bins(),
            reference_sequence.metadata(),
        )?;
    }

    Ok(())
}

fn write_bins<W>(
    writer: &mut W,
    depth: u8,
    bins: &[Bin],
    metadata: Option<&Metadata>,
) -> io::Result<()>
where
    W: Write,
{
    let mut n_bin =
        i32::try_from(bins.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if metadata.is_some() {
        n_bin += 1;
    }

    writer.write_i32::<LittleEndian>(n_bin)?;

    for bin in bins {
        let bin_id =
            u32::try_from(bin.id()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        writer.write_u32::<LittleEndian>(bin_id)?;

        let loffset = u64::from(bin.loffset());
        writer.write_u64::<LittleEndian>(loffset)?;

        write_chunks(writer, bin.chunks())?;
    }

    if let Some(m) = metadata {
        write_metadata(writer, depth, m)?;
    }

    Ok(())
}

fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: Write,
{
    let n_chunk =
        i32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_chunk)?;

    for chunk in chunks {
        let chunk_beg = u64::from(chunk.start());
        writer.write_u64::<LittleEndian>(chunk_beg)?;

        let chunk_end = u64::from(chunk.start());
        writer.write_u64::<LittleEndian>(chunk_end)?;
    }

    Ok(())
}

fn write_metadata<W>(writer: &mut W, depth: u8, metadata: &Metadata) -> io::Result<()>
where
    W: Write,
{
    const N_CHUNK: i32 = 2;

    let bin_id = u32::try_from(Bin::metadata_id(depth))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(bin_id)?;

    let loffset = u64::from(bgzf::VirtualPosition::default());
    writer.write_u64::<LittleEndian>(loffset)?;

    writer.write_i32::<LittleEndian>(N_CHUNK)?;

    let ref_beg = u64::from(metadata.start_position());
    writer.write_u64::<LittleEndian>(ref_beg)?;

    let ref_end = u64::from(metadata.end_position());
    writer.write_u64::<LittleEndian>(ref_end)?;

    let n_mapped = metadata.mapped_record_count();
    writer.write_u64::<LittleEndian>(n_mapped)?;

    let n_unmapped = metadata.unmapped_record_count();
    writer.write_u64::<LittleEndian>(n_unmapped)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_metadata() -> io::Result<()> {
        let mut buf = Vec::new();
        let depth = 5;
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        write_metadata(&mut buf, depth, &metadata)?;

        let expected = [
            0x4a, 0x92, 0x00, 0x00, // bin = 37450
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // loffset = 0
            0x02, 0x00, 0x00, 0x00, // chunks = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
