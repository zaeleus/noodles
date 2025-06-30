use std::{
    io::{self, Read},
    mem,
};

use crate::gzi::Index;

pub(super) fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: Read,
{
    let len = read_u64_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let offsets: Vec<_> = (0..len)
        .map(|_| Ok((read_u64_le(reader)?, read_u64_le(reader)?)))
        .collect::<io::Result<_>>()?;

    match read_u8(reader) {
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected trailing data",
        )),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(Index::from(offsets)),
        Err(e) => Err(e),
    }
}

fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u8>()];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u64_le<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u64>()];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_index() -> io::Result<()> {
        let src = [
            0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // len = 2
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x02, 0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 23810
            0x01, 0x52, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 86529
        ];

        let mut reader = &src[..];
        assert_eq!(
            read_index(&mut reader)?,
            Index::from(vec![(4668, 21294), (23810, 86529)])
        );

        Ok(())
    }

    #[test]
    fn test_read_index_with_no_entries() -> io::Result<()> {
        let src = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // len = 0
        let mut reader = &src[..];
        assert_eq!(read_index(&mut reader)?, Index::default());
        Ok(())
    }

    #[test]
    fn test_read_index_with_fewer_than_len_entries() -> io::Result<()> {
        let src = [0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // len = 1
        let mut reader = &src[..];

        assert!(matches!(
            read_index(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    #[test]
    fn test_read_index_with_trailing_data() -> io::Result<()> {
        let src = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // len = 0
            0x00,
        ];

        let mut reader = &src[..];

        assert!(matches!(
            read_index(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
