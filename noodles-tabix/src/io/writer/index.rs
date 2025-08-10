mod header;
mod magic_number;
mod reference_sequences;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_csi::BinningIndex;

use self::{
    header::write_header, magic_number::write_magic_number,
    reference_sequences::write_reference_sequences,
};
use super::num::write_i32_le;
use crate::Index;

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;

    let reference_sequence_count = i32::try_from(index.reference_sequences().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, reference_sequence_count)?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tabix header"))?;
    write_header(writer, header)?;

    write_reference_sequences(writer, index.reference_sequences())?;

    if let Some(n) = index.unplaced_unmapped_record_count() {
        writer.write_u64::<LittleEndian>(n)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;
    use noodles_bgzf as bgzf;
    use noodles_csi::binning_index::index::{
        Header, ReferenceSequence,
        reference_sequence::{Bin, bin::Chunk},
    };

    use super::*;
    use crate::{MAGIC_NUMBER, io::writer::num::write_u32_le};

    #[test]
    fn test_write_index() -> io::Result<()> {
        let chunks = vec![Chunk::new(
            bgzf::VirtualPosition::from(509268599425),
            bgzf::VirtualPosition::from(509268599570),
        )];
        let bins = [(16385, Bin::new(chunks))].into_iter().collect();
        let intervals = vec![bgzf::VirtualPosition::from(337)];
        let references = vec![ReferenceSequence::new(bins, intervals, None)];

        let reference_sequence_names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();

        let header = Header::builder()
            .set_reference_sequence_names(reference_sequence_names)
            .build();

        let index = Index::builder()
            .set_header(header)
            .set_reference_sequences(references)
            .build();

        let mut buf = Vec::new();
        write_index(&mut buf, &index)?;

        let mut expected = Vec::new();
        expected.write_all(&MAGIC_NUMBER)?; // magic
        write_i32_le(&mut expected, 1)?; // n_ref
        write_i32_le(&mut expected, 0)?; // format
        write_i32_le(&mut expected, 1)?; // col_seq
        write_i32_le(&mut expected, 4)?; // col_beg
        write_i32_le(&mut expected, 5)?; // col_end
        write_i32_le(&mut expected, i32::from(b'#'))?; // meta
        write_i32_le(&mut expected, 0)?; // skip
        write_i32_le(&mut expected, 8)?; // l_nm
        expected.write_all(b"sq0\x00sq1\x00")?; // names
        write_u32_le(&mut expected, 1)?; // n_bin
        write_u32_le(&mut expected, 16385)?; // bin
        write_u32_le(&mut expected, 1)?; // n_chunk
        expected.write_u64::<LittleEndian>(509268599425)?; // chunk_beg
        expected.write_u64::<LittleEndian>(509268599570)?; // chunk_end
        write_u32_le(&mut expected, 1)?; // n_intv
        expected.write_u64::<LittleEndian>(337)?; // ioffset

        assert_eq!(buf, expected);

        Ok(())
    }
}
