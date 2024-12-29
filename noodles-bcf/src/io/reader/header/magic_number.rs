use std::io::{self, Read};

pub(crate) fn read_magic_number<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    use crate::MAGIC_NUMBER;

    let mut buf = [0; 3];
    reader.read_exact(&mut buf)?;

    if buf == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF header",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic_number() {
        let data = b"BCF";
        let mut reader = &data[..];
        assert!(read_magic_number(&mut reader).is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAM";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
