use std::io;

use crate::io::reader::num::{read_u8, read_u32_le};

#[derive(Debug, Eq, PartialEq)]
pub(super) enum CompressionMethod {
    RansNx16,
    AdaptiveArithmeticCoder,
}

pub(super) fn read_header(src: &mut &[u8]) -> io::Result<(usize, usize, CompressionMethod)> {
    let ulen = read_u32_le(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_names = read_u32_le(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let compression_method = read_compression_method(src)?;

    Ok((ulen, n_names, compression_method))
}

fn read_compression_method(src: &mut &[u8]) -> io::Result<CompressionMethod> {
    read_u8(src).map(|n| match n {
        0 => CompressionMethod::RansNx16,
        _ => CompressionMethod::AdaptiveArithmeticCoder,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_compression_method() -> io::Result<()> {
        fn t(mut src: &[u8], expected: CompressionMethod) -> io::Result<()> {
            assert_eq!(read_compression_method(&mut src)?, expected);
            Ok(())
        }

        t(&[0x00], CompressionMethod::RansNx16)?;
        t(&[0x01], CompressionMethod::AdaptiveArithmeticCoder)?;
        t(&[0x02], CompressionMethod::AdaptiveArithmeticCoder)?;
        t(&[0xff], CompressionMethod::AdaptiveArithmeticCoder)?;

        Ok(())
    }
}
