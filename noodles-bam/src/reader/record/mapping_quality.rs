use std::{io, mem};

use bytes::Buf;
use noodles_sam::{self as sam, record::MappingQuality};

pub fn get_mapping_quality<B>(src: &mut B) -> io::Result<Option<MappingQuality>>
where
    B: Buf,
{
    use sam::record::mapping_quality::MISSING;

    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_u8() {
        MISSING => Ok(None),
        n => Ok(MappingQuality::new(n)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_mapping_quality() -> io::Result<()> {
        fn t(mut buf: &[u8], expected: Option<MappingQuality>) -> io::Result<()> {
            let actual = get_mapping_quality(&mut buf)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0xff], None)?;
        t(&[0x08], MappingQuality::new(8))?;

        Ok(())
    }
}
