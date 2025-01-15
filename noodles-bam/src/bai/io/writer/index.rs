use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_csi::{
    binning_index::{
        index::{
            reference_sequence::{bin::Chunk, index::LinearIndex, Bin, Metadata},
            ReferenceSequence,
        },
        ReferenceSequence as _,
    },
    BinningIndex,
};

use crate::bai::{Index, MAGIC_NUMBER};

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&MAGIC_NUMBER)?;

    let n_ref = u32::try_from(index.reference_sequences().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_ref)?;

    for reference_sequence in index.reference_sequences() {
        write_reference_sequence(writer, reference_sequence)?;
    }

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        writer.write_u64::<LittleEndian>(n_no_coor)?;
    }

    Ok(())
}

fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &ReferenceSequence<LinearIndex>,
) -> io::Result<()>
where
    W: Write,
{
    let mut n_bin = u32::try_from(reference_sequence.bins().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if reference_sequence.metadata().is_some() {
        n_bin = n_bin
            .checked_add(1)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))?;
    }

    writer.write_u32::<LittleEndian>(n_bin)?;

    for (&id, bin) in reference_sequence.bins() {
        write_bin(writer, id, bin)?;
    }

    if let Some(metadata) = reference_sequence.metadata() {
        write_metadata(writer, metadata)?;
    }

    let n_intv = u32::try_from(reference_sequence.index().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_intv)?;

    for interval in reference_sequence.index() {
        let ioffset = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioffset)?;
    }

    Ok(())
}

fn write_bin<W>(writer: &mut W, id: usize, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(id)?;

    let n_chunk = u32::try_from(bin.chunks().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_chunk)?;

    for chunk in bin.chunks() {
        write_chunk(writer, chunk)?;
    }

    Ok(())
}

fn write_chunk<W>(writer: &mut W, chunk: &Chunk) -> io::Result<()>
where
    W: Write,
{
    let chunk_beg = u64::from(chunk.start());
    writer.write_u64::<LittleEndian>(chunk_beg)?;

    let chunk_end = u64::from(chunk.end());
    writer.write_u64::<LittleEndian>(chunk_end)?;

    Ok(())
}

fn write_metadata<W>(writer: &mut W, metadata: &Metadata) -> io::Result<()>
where
    W: Write,
{
    use crate::bai::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);
    const METADATA_CHUNK_COUNT: usize = 2;

    let id =
        u32::try_from(METADATA_ID).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(id)?;

    let n_chunk = u32::try_from(METADATA_CHUNK_COUNT)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_chunk)?;

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
    use noodles_bgzf as bgzf;

    use super::*;

    #[test]
    fn test_write_index() -> io::Result<()> {
        let chunks = vec![Chunk::new(
            bgzf::VirtualPosition::from(509268599425),
            bgzf::VirtualPosition::from(509268599570),
        )];
        let bins = [(16385, Bin::new(chunks))].into_iter().collect();
        let intervals = vec![bgzf::VirtualPosition::from(337)];
        let reference_sequences = vec![ReferenceSequence::new(bins, intervals, None)];
        let index = Index::builder()
            .set_reference_sequences(reference_sequences)
            .build();

        let mut buf = Vec::new();
        write_index(&mut buf, &index)?;

        let mut expected = Vec::new();
        // magic
        expected.write_all(&MAGIC_NUMBER)?;
        // n_ref
        expected.write_u32::<LittleEndian>(1)?;
        // n_bin
        expected.write_u32::<LittleEndian>(1)?;
        // bin
        expected.write_u32::<LittleEndian>(16385)?;
        // n_chunk
        expected.write_u32::<LittleEndian>(1)?;
        // chunk_beg
        expected.write_u64::<LittleEndian>(509268599425)?;
        // chunk_end
        expected.write_u64::<LittleEndian>(509268599570)?;
        // n_intv
        expected.write_u32::<LittleEndian>(1)?;
        // ioffset
        expected.write_u64::<LittleEndian>(337)?;

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_metadata() -> io::Result<()> {
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        let mut buf = Vec::new();
        write_metadata(&mut buf, &metadata)?;

        let expected = [
            0x4a, 0x92, 0x00, 0x00, // bin = 37450
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
