use std::{io, num::NonZero};

use crate::io::reader::num::{read_u8, read_u32_le};

#[derive(Debug, Eq, PartialEq)]
pub(super) enum CompressionMethod {
    RansNx16,
    AdaptiveArithmeticCoder,
}

pub(super) fn read_header(
    src: &mut &[u8],
) -> io::Result<(usize, NonZero<usize>, CompressionMethod)> {
    let uncompressed_size = read_uncompressed_size(src)?;
    let name_count = read_name_count(src)?;
    let compression_method = read_compression_method(src)?;
    Ok((uncompressed_size, name_count, compression_method))
}

fn read_uncompressed_size(src: &mut &[u8]) -> io::Result<usize> {
    read_u32_le(src)
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

fn read_name_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_u32_le(src).and_then(|n| {
        usize::try_from(n)
            .and_then(NonZero::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_compression_method(src: &mut &[u8]) -> io::Result<CompressionMethod> {
    read_u8(src).and_then(decode_compression_method)
}

// § "Name tokenisation codec" (2025-04-07): "`use_arith`: Whether compression is arithmetic (1) or
// rANS Nx16 (0)".
fn decode_compression_method(n: u8) -> io::Result<CompressionMethod> {
    match n {
        0 => Ok(CompressionMethod::RansNx16),
        1 => Ok(CompressionMethod::AdaptiveArithmeticCoder),
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
    fn test_read_name_count() -> io::Result<()> {
        assert_eq!(
            read_name_count(&mut &[0x01, 0x00, 0x00, 0x00][..])?,
            NonZero::<usize>::MIN
        );

        assert!(matches!(
            read_name_count(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_name_count(&mut &[0x00, 0x00, 0x00, 0x00][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_compression_method() -> io::Result<()> {
        fn t(mut src: &[u8], expected: CompressionMethod) -> io::Result<()> {
            assert_eq!(read_compression_method(&mut src)?, expected);
            Ok(())
        }

        t(&[0x00], CompressionMethod::RansNx16)?;
        t(&[0x01], CompressionMethod::AdaptiveArithmeticCoder)?;

        assert!(matches!(
            read_compression_method(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_compression_method(&mut &[0x02][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
