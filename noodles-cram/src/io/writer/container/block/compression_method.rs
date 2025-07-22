use std::io::{self, Write};

use crate::{container::block::CompressionMethod, io::writer::num::write_u8};

pub(super) fn write_compression_method<W>(
    writer: &mut W,
    compression_method: CompressionMethod,
) -> io::Result<()>
where
    W: Write,
{
    let n = encode(compression_method);
    write_u8(writer, n)
}

fn encode(compression_method: CompressionMethod) -> u8 {
    match compression_method {
        CompressionMethod::None => 0,
        CompressionMethod::Gzip => 1,
        CompressionMethod::Bzip2 => 2,
        CompressionMethod::Lzma => 3,
        CompressionMethod::Rans4x8 => 4,
        CompressionMethod::RansNx16 => 5,
        CompressionMethod::AdaptiveArithmeticCoding => 6,
        CompressionMethod::Fqzcomp => 7,
        CompressionMethod::NameTokenizer => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_compression_method() -> io::Result<()> {
        let mut buf = Vec::new();
        write_compression_method(&mut buf, CompressionMethod::None)?;
        assert_eq!(buf, [0x00]);
        Ok(())
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(CompressionMethod::None), 0);
        assert_eq!(encode(CompressionMethod::Gzip), 1);
        assert_eq!(encode(CompressionMethod::Bzip2), 2);
        assert_eq!(encode(CompressionMethod::Lzma), 3);
        assert_eq!(encode(CompressionMethod::Rans4x8), 4);
        assert_eq!(encode(CompressionMethod::RansNx16), 5);
        assert_eq!(encode(CompressionMethod::AdaptiveArithmeticCoding), 6);
        assert_eq!(encode(CompressionMethod::Fqzcomp), 7);
        assert_eq!(encode(CompressionMethod::NameTokenizer), 8);
    }
}
