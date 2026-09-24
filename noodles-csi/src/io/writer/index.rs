pub(crate) mod header;
mod magic_number;
mod reference_sequences;

use std::io::{self, Write};

pub use self::header::write_header;
use self::{
    header::write_aux, magic_number::write_magic_number,
    reference_sequences::write_reference_sequences,
};
use super::num::{write_i32_le, write_u64_le};
use crate::{BinningIndex, Index};

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;

    let min_shift = i32::from(index.min_shift());
    write_i32_le(writer, min_shift)?;

    write_depth(writer, index.depth())?;

    write_aux(writer, index.header())?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences())?;

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        write_u64_le(writer, n_no_coor)?;
    }

    Ok(())
}

fn write_depth<W>(writer: &mut W, depth: u8) -> io::Result<()>
where
    W: Write,
{
    let n = i32::from(depth);
    write_i32_le(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_depth() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, depth: u8, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_depth(buf, depth)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00, 0x00, 0x00, 0x00])?;
        t(&mut buf, 5, &[0x05, 0x00, 0x00, 0x00])?;
        t(&mut buf, 9, &[0x09, 0x00, 0x00, 0x00])?;

        t(&mut buf, 10, &[0x0a, 0x00, 0x00, 0x00])?; // FIXME

        Ok(())
    }
}
