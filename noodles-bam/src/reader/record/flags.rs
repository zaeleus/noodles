use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::Flags;

pub(crate) fn get_flags<B>(src: &mut B) -> io::Result<Flags>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(Flags::from(src.get_u16_le()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_flags() -> io::Result<()> {
        let mut src = &[0x00, 0x00][..];
        assert_eq!(get_flags(&mut src)?, Flags::empty());

        let mut src = &[0x04, 0x00][..];
        assert_eq!(get_flags(&mut src)?, Flags::UNMAPPED);

        let mut src = &[][..];
        assert!(matches!(
            get_flags(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &[0x00][..];
        assert!(matches!(
            get_flags(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
