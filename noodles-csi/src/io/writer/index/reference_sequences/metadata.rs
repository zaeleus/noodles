use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;

use crate::{
    binning_index::index::reference_sequence::{Bin, Metadata},
    io::writer::num::write_i32_le,
};

pub(super) fn write_metadata<W>(writer: &mut W, depth: u8, metadata: &Metadata) -> io::Result<()>
where
    W: Write,
{
    const N_CHUNK: i32 = 2;

    let bin_id = u32::try_from(Bin::metadata_id(depth))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(bin_id)?;

    let loffset = u64::from(bgzf::VirtualPosition::default());
    writer.write_u64::<LittleEndian>(loffset)?;

    write_i32_le(writer, N_CHUNK)?;

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
