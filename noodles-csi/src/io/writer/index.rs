pub(crate) mod header;
mod reference_sequences;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use self::{header::write_aux, reference_sequences::write_reference_sequences};
use super::num::write_i32_le;
use crate::{BinningIndex, Index, io::MAGIC_NUMBER};

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic(writer)?;

    let min_shift = i32::from(index.min_shift());
    write_i32_le(writer, min_shift)?;

    let depth = i32::from(index.depth());
    write_i32_le(writer, depth)?;

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
    writer.write_all(&MAGIC_NUMBER)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_magic() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic(&mut buf)?;
        assert_eq!(buf, MAGIC_NUMBER);
        Ok(())
    }
}
