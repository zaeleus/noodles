use std::io;

use byteorder::{LittleEndian, ReadBytesExt};

pub(crate) fn read_qual(src: &mut &[u8]) -> io::Result<Option<f32>> {
    use crate::record::codec::value::Float;

    match src.read_f32::<LittleEndian>().map(Float::from)? {
        Float::Value(n) => Ok(Some(n)),
        Float::Missing => Ok(None),
        qual => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid qual: {qual:?}"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_qual() -> io::Result<()> {
        let mut src = &[0x01, 0x00, 0x80, 0x7f][..];
        assert!(read_qual(&mut src)?.is_none());

        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_qual(&mut src)?, Some(0.0));

        let mut src = &[0x02, 0x00, 0x80, 0x7f][..];
        assert!(matches!(
            read_qual(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
