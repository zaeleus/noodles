use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::gzi::Index;

pub(super) fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    let index = index.as_ref();

    let len =
        u64::try_from(index.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u64::<LittleEndian>(len)?;

    for &(compressed_pos, uncompressed_pos) in index {
        writer.write_u64::<LittleEndian>(compressed_pos)?;
        writer.write_u64::<LittleEndian>(uncompressed_pos)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_index() -> io::Result<()> {
        let mut buf = Vec::new();

        let index = Index::from(vec![(4668, 21294), (23810, 86529)]);
        write_index(&mut buf, &index)?;

        let expected = [
            0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // len = 2
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x02, 0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 23810
            0x01, 0x52, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 86529
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_index_with_no_entries() -> io::Result<()> {
        let mut buf = Vec::new();
        let index = Index::default();
        write_index(&mut buf, &index)?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);
        Ok(())
    }
}
