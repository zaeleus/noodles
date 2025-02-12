use std::io;

use bytes::Buf;

use crate::container::block::CompressionMethod;

pub(super) fn get_compression_method<B>(src: &mut B) -> io::Result<CompressionMethod>
where
    B: Buf,
{
    let n = src
        .try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))?;

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
        fn t(mut src: &[u8], expected: CompressionMethod) -> io::Result<()> {
            let actual = get_compression_method(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], CompressionMethod::None)?;
        t(&[0x01], CompressionMethod::Gzip)?;
        t(&[0x02], CompressionMethod::Bzip2)?;
        t(&[0x03], CompressionMethod::Lzma)?;
        t(&[0x04], CompressionMethod::Rans4x8)?;
        t(&[0x05], CompressionMethod::RansNx16)?;
        t(&[0x06], CompressionMethod::AdaptiveArithmeticCoding)?;
        t(&[0x07], CompressionMethod::Fqzcomp)?;
        t(&[0x08], CompressionMethod::NameTokenizer)?;

        let mut src = &[][..];
        assert!(matches!(
            get_compression_method(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &[0x09][..];
        assert!(matches!(
            get_compression_method(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
