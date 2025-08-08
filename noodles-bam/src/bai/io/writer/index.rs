mod magic_number;
mod reference_sequences;

use std::io::{self, Write};

use noodles_csi::BinningIndex;

use self::{magic_number::write_magic_number, reference_sequences::write_reference_sequences};
use crate::{bai::Index, io::writer::num::write_u64_le};

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;
    write_reference_sequences(writer, index.reference_sequences())?;

    if let Some(n) = index.unplaced_unmapped_record_count() {
        write_unplaced_unmapped_record_count(writer, n)?;
    }

    Ok(())
}

fn write_unplaced_unmapped_record_count<W>(
    writer: &mut W,
    unplaced_unmapped_record_count: u64,
) -> io::Result<()>
where
    W: Write,
{
    write_u64_le(writer, unplaced_unmapped_record_count)
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;
    use noodles_csi::binning_index::index::{
        ReferenceSequence,
        reference_sequence::{Bin, bin::Chunk},
    };

    use super::*;
    use crate::{bai::MAGIC_NUMBER, io::writer::num::write_u32_le};

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
        expected.write_all(&MAGIC_NUMBER)?; // magic
        write_u32_le(&mut expected, 1)?; // n_ref
        write_u32_le(&mut expected, 1)?; // n_bin
        write_u32_le(&mut expected, 16385)?; // bin
        write_u32_le(&mut expected, 1)?; // n_chunk
        write_u64_le(&mut expected, 509268599425)?; // chunk_beg
        write_u64_le(&mut expected, 509268599570)?; // chunk_end
        write_u32_le(&mut expected, 1)?; // n_intv
        write_u64_le(&mut expected, 337)?; // ioffset

        assert_eq!(buf, expected);

        Ok(())
    }
}
