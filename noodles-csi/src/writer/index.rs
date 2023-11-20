pub(crate) mod header;
mod reference_sequences;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use self::{header::write_aux, reference_sequences::write_reference_sequences};
use crate::{index::reference_sequence::index::BinnedIndex, BinningIndex, Index};

pub(super) fn write_index<W>(writer: &mut W, index: &Index<BinnedIndex>) -> io::Result<()>
where
    W: Write,
{
    write_magic(writer)?;

    let min_shift = i32::from(index.min_shift());
    writer.write_i32::<LittleEndian>(min_shift)?;

    let depth = i32::from(index.depth());
    writer.write_i32::<LittleEndian>(depth)?;

    write_aux(writer, index.header())?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences())?;

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        writer.write_u64::<LittleEndian>(n_no_coor)?;
    }

    Ok(())
}

fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    use crate::io::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_magic() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic(&mut buf)?;
        let expected = b"CSI\x01";
        assert_eq!(buf, expected);
        Ok(())
    }
}
