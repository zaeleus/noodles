use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::record::data::field::value::Subtype;

pub fn read_subtype<R>(reader: &mut R) -> io::Result<Subtype>
where
    R: Read,
{
    reader.read_u8().and_then(|n| {
        Subtype::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_subtype() -> io::Result<()> {
        let data = [b'i'];
        let mut reader = &data[..];
        assert_eq!(read_subtype(&mut reader)?, Subtype::Int32);

        let data = [b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            read_subtype(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
