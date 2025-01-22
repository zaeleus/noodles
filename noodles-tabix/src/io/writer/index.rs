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
use crate::Index;

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;

    let reference_sequence_count = i32::try_from(index.reference_sequences().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(reference_sequence_count)?;

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
        reference_sequence::{bin::Chunk, Bin},
        Header, ReferenceSequence,
    };

    use super::*;
    use crate::MAGIC_NUMBER;

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
        expected.write_i32::<LittleEndian>(1)?; // n_ref
        expected.write_i32::<LittleEndian>(0)?; // format
        expected.write_i32::<LittleEndian>(1)?; // col_seq
        expected.write_i32::<LittleEndian>(4)?; // col_beg
        expected.write_i32::<LittleEndian>(5)?; // col_end
        expected.write_i32::<LittleEndian>(i32::from(b'#'))?; // meta
        expected.write_i32::<LittleEndian>(0)?; // skip
        expected.write_i32::<LittleEndian>(8)?; // l_nm
        expected.write_all(b"sq0\x00sq1\x00")?; // names
        expected.write_u32::<LittleEndian>(1)?; // n_bin
        expected.write_u32::<LittleEndian>(16385)?; // bin
        expected.write_u32::<LittleEndian>(1)?; // n_chunk
        expected.write_u64::<LittleEndian>(509268599425)?; // chunk_beg
        expected.write_u64::<LittleEndian>(509268599570)?; // chunk_end
        expected.write_u32::<LittleEndian>(1)?; // n_intv
        expected.write_u64::<LittleEndian>(337)?; // ioffset

        assert_eq!(buf, expected);

        Ok(())
    }
}
