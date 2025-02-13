use std::io;

use bytes::Buf;

use crate::container::block::CompressionMethod;

pub(super) fn get_compression_method<B>(src: &mut B) -> io::Result<CompressionMethod>
where
    B: Buf,
{
    src.try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))
        .and_then(decode)
}

fn decode(n: u8) -> io::Result<CompressionMethod> {
    match n {
        0 => Ok(CompressionMethod::None),
        1 => Ok(CompressionMethod::Gzip),
        2 => Ok(CompressionMethod::Bzip2),
        3 => Ok(CompressionMethod::Lzma),
        4 => Ok(CompressionMethod::Rans4x8),
        5 => Ok(CompressionMethod::RansNx16),
        6 => Ok(CompressionMethod::AdaptiveArithmeticCoding),
        7 => Ok(CompressionMethod::Fqzcomp),
        8 => Ok(CompressionMethod::NameTokenizer),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid compression method",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_compression_method() -> io::Result<()> {
        assert_eq!(
            get_compression_method(&mut &[0x00][..])?,
            CompressionMethod::None
        );

        assert!(matches!(
            get_compression_method(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        assert_eq!(decode(0)?, CompressionMethod::None);
        assert_eq!(decode(1)?, CompressionMethod::Gzip);
        assert_eq!(decode(2)?, CompressionMethod::Bzip2);
        assert_eq!(decode(3)?, CompressionMethod::Lzma);
        assert_eq!(decode(4)?, CompressionMethod::Rans4x8);
        assert_eq!(decode(5)?, CompressionMethod::RansNx16);
        assert_eq!(decode(6)?, CompressionMethod::AdaptiveArithmeticCoding);
        assert_eq!(decode(7)?, CompressionMethod::Fqzcomp);
        assert_eq!(decode(8)?, CompressionMethod::NameTokenizer);

        assert!(matches!(
            decode(9),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
