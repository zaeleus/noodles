use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::container::block::CompressionMethod;

pub(super) fn write_compression_method<W>(
    writer: &mut W,
    compression_method: CompressionMethod,
) -> io::Result<()>
where
    W: Write,
{
    let n = match compression_method {
        CompressionMethod::None => 0,
        CompressionMethod::Gzip => 1,
        CompressionMethod::Bzip2 => 2,
        CompressionMethod::Lzma => 3,
        CompressionMethod::Rans4x8 => 4,
        CompressionMethod::RansNx16 => 5,
        CompressionMethod::AdaptiveArithmeticCoding => 6,
        CompressionMethod::Fqzcomp => 7,
        CompressionMethod::NameTokenizer => 8,
    };

    writer.write_u8(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_compression_method() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            compression_method: CompressionMethod,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_compression_method(buf, compression_method)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, CompressionMethod::None, &[0x00])?;
        t(&mut buf, CompressionMethod::Gzip, &[0x01])?;
        t(&mut buf, CompressionMethod::Bzip2, &[0x02])?;
        t(&mut buf, CompressionMethod::Lzma, &[0x03])?;
        t(&mut buf, CompressionMethod::Rans4x8, &[0x04])?;
        t(&mut buf, CompressionMethod::RansNx16, &[0x05])?;
        t(
            &mut buf,
            CompressionMethod::AdaptiveArithmeticCoding,
            &[0x06],
        )?;
        t(&mut buf, CompressionMethod::Fqzcomp, &[0x07])?;
        t(&mut buf, CompressionMethod::NameTokenizer, &[0x08])?;

        Ok(())
    }
}
