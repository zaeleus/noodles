use std::io::{self, Read};

use crate::io::MAGIC_NUMBER;

type Buf = [u8; MAGIC_NUMBER.len()];

pub(super) fn read_magic_number<R>(reader: &mut R) -> io::Result<Buf>
where
    R: Read,
{
    let mut buf = [0; MAGIC_NUMBER.len()];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

pub(crate) fn validate(buf: Buf) -> io::Result<()> {
    if is_valid(buf) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

fn is_valid(buf: Buf) -> bool {
    buf == MAGIC_NUMBER
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic_number() -> io::Result<()> {
        let src = b"BAM\x01";
        let mut reader = &src[..];
        assert_eq!(read_magic_number(&mut reader)?, *b"BAM\x01");

        let mut src = &[][..];
        assert!(matches!(
            read_magic_number(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(*b"BAM\x01"));
        assert!(!is_valid([0x00, 0x00, 0x00, 0x00]));
    }
}
