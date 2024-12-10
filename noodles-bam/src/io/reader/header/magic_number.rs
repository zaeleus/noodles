use std::io::{self, Read};

use super::Reader;
use crate::MAGIC_NUMBER;

pub(super) fn read_magic_number<R>(reader: &mut Reader<R>) -> io::Result<()>
where
    R: Read,
{
    let magic_number = reader.read_magic_number()?;

    if magic_number == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic_number() -> io::Result<()> {
        let mut src = &b"BAM\x01"[..];
        let mut reader = Reader::new(&mut src);
        assert!(read_magic_number(&mut reader).is_ok());

        let mut src = &[][..];
        let mut reader = Reader::new(&mut src);
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"MThd"[..];
        let mut reader = Reader::new(&mut src);
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
