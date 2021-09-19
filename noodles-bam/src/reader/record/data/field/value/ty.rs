use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::ReadBytesExt;

use crate::record::data::field::value::Type;

pub fn read_type<R>(reader: &mut R) -> io::Result<Type>
where
    R: Read,
{
    reader
        .read_u8()
        .and_then(|n| Type::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() -> io::Result<()> {
        let data = [b'i'];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int32);

        let data = [b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            read_type(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
