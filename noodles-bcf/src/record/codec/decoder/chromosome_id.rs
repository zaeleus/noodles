use std::io;

use byteorder::{LittleEndian, ReadBytesExt};

pub(crate) fn read_chrom(src: &mut &[u8]) -> io::Result<usize> {
    src.read_i32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_chrom() -> io::Result<()> {
        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_chrom(&mut src)?, 0);

        let mut src = &[][..];
        assert!(matches!(
            read_chrom(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &[0xff, 0xff, 0xff, 0xff][..];
        assert!(matches!(
            read_chrom(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
