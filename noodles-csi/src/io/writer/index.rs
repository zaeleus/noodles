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
use crate::{BinningIndex, Index, io::MAX_DEPTH};

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;

    write_min_shift(writer, index.min_shift())?;
    write_depth(writer, index.depth())?;

    write_aux(writer, index.header())?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences())?;

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        write_u64_le(writer, n_no_coor)?;
    }

    Ok(())
}

fn write_min_shift<W>(writer: &mut W, min_shift: u8) -> io::Result<()>
where
    W: Write,
{
    let n = i32::from(min_shift);
    write_i32_le(writer, n)
}

fn write_depth<W>(writer: &mut W, depth: u8) -> io::Result<()>
where
    W: Write,
{
    if depth > MAX_DEPTH {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid depth"));
    }

    let n = i32::from(depth);
    write_i32_le(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_min_shift() -> io::Result<()> {
        let mut buf = Vec::new();
        write_min_shift(&mut buf, 14)?;
        assert_eq!(buf, [0x0e, 0x00, 0x00, 0x00]);
        Ok(())
    }

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

        buf.clear();
        assert!(matches!(
            write_depth(&mut buf, 10),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
