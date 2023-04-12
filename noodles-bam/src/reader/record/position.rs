use std::{io, mem};

use bytes::Buf;
use noodles_core::Position;

pub(crate) fn get_position<B>(src: &mut B) -> io::Result<Option<Position>>
where
    B: Buf,
{
    const MISSING: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_i32_le() {
        MISSING => Ok(None),
        n => usize::try_from(n + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Position::new),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_position() -> io::Result<()> {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert!(get_position(&mut src)?.is_none());

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_position(&mut src)?, Some(Position::MIN));

        let mut src = &[][..];
        assert!(matches!(
            get_position(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            get_position(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
