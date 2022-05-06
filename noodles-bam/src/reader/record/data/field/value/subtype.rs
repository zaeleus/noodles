use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::data::field::value::Subtype;

pub fn get_subtype<B>(src: &mut B) -> io::Result<Subtype>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Subtype::try_from(src.get_u8()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_subtype() -> io::Result<()> {
        let data = [b'i'];
        let mut reader = &data[..];
        assert_eq!(get_subtype(&mut reader)?, Subtype::Int32);

        let data = [b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            get_subtype(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
