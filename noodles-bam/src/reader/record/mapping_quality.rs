use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::MappingQuality;

pub fn get_mapping_quality<B>(src: &mut B) -> io::Result<Option<MappingQuality>>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(MappingQuality::new(src.get_u8()))
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
