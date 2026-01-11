mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::{
    borrow::Cow,
    io::{self, Write},
    num::NonZero,
};

use super::Flags;
use crate::io::writer::num::{write_u8, write_uint7};

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    use crate::codecs::rans_nx16::encode::bit_pack;

    let mut src = Cow::from(src);
    let mut dst = Vec::new();

    write_flags(&mut dst, flags)?;

    if flags.has_uncompressed_size() {
        write_uncompressed_size(&mut dst, src.len())?;
    }

    if flags.is_striped() {
        let buf = stripe::encode(&src)?;
        dst.extend(buf);
        return Ok(dst);
    }

    if flags.is_bit_packed() {
        match bit_pack::build_context(&src) {
            Ok(ctx) => {
                src = Cow::from(bit_pack::encode(&src, &ctx));
                bit_pack::write_context(&mut dst, &ctx, src.len())?;
            }
            Err(
                bit_pack::context::BuildContextError::EmptyAlphabet
                | bit_pack::context::BuildContextError::TooManySymbols(_),
            ) => {
                flags.remove(Flags::PACK);
                dst[0] = u8::from(flags);
            }
        }
    }

    if flags.is_uncompressed() {
        dst.write_all(&src)?;
    } else if flags.uses_external_codec() {
        encode_ext(&src, &mut dst)?;
    } else if flags.is_rle() {
        rle::encode(&src, flags, &mut dst)?;
    } else if flags.order() == 0 {
        order_0::encode(&src, &mut dst)?;
    } else {
        order_1::encode(&src, &mut dst)?;
    }

    Ok(dst)
}

fn write_flags(dst: &mut Vec<u8>, flags: Flags) -> io::Result<()> {
    write_u8(dst, u8::from(flags))
}

fn write_uncompressed_size(dst: &mut Vec<u8>, uncompressed_size: usize) -> io::Result<()> {
    let n = u32::try_from(uncompressed_size)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_uint7(dst, n)
}

fn write_symbol_count(dst: &mut Vec<u8>, symbol_count: NonZero<usize>) -> io::Result<()> {
    let n = u8::try_from(usize::from(symbol_count))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_u8(dst, n)
}

fn encode_ext(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    use bzip2::write::BzEncoder;

    let mut encoder = BzEncoder::new(dst, Default::default());
    encoder.write_all(src)?;
    encoder.finish()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_ext() -> io::Result<()> {
        use crate::codecs::bzip2;

        let actual = encode(Flags::EXT, b"noodles")?;

        let mut expected = vec![0x04, 0x07];

        let compression_level = ::bzip2::Compression::default();
        let data = bzip2::encode(compression_level, b"noodles")?;
        expected.extend(data);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_stripe() -> io::Result<()> {
        let actual = encode(Flags::STRIPE, b"noodles")?;

        let expected = [
            0x08, 0x07, 0x04, 0x09, 0x09, 0x09, 0x08, 0x00, 0x02, 0x6f, 0x00, 0xff, 0xa7, 0xab,
            0x62, 0x00, 0x00, 0x02, 0x70, 0x00, 0xff, 0x84, 0x92, 0x1b, 0x00, 0x00, 0x02, 0x74,
            0x00, 0xf7, 0x27, 0xdb, 0x24, 0x00, 0x00, 0x01, 0x65, 0x00, 0xfd, 0x77, 0x20, 0xb0,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_0() -> io::Result<()> {
        let actual = encode(Flags::empty(), b"noodles")?;

        let expected = [
            0x00, 0x07, 0x74, 0x00, 0xf4, 0xe5, 0xb7, 0x4e, 0x50, 0x0f, 0x2e, 0x97, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER, b"noodles")?;

        let expected = [
            0x01, 0x07, 0x74, 0x00, 0xf4, 0xe3, 0x83, 0x41, 0xe2, 0x9a, 0xef, 0x53, 0x50, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_cat() -> io::Result<()> {
        let actual = encode(Flags::CAT, b"noodles")?;
        let expected = [0x20, 0x07, 0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73];
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_encode_rle_with_order_0() -> io::Result<()> {
        let actual = encode(Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x40, 0x0d, 0x74, 0x00, 0xf3, 0x4b, 0x21, 0x10, 0xa8, 0xe3, 0x84, 0xfe, 0x6b, 0x22,
            0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_rle_with_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER | Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x41, 0x0d, 0x74, 0x00, 0xf3, 0x4a, 0x89, 0x79, 0xc1, 0xe8, 0xc3, 0xc5, 0x62, 0x31,
            0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_pack() -> io::Result<()> {
        let actual = encode(Flags::CAT | Flags::PACK, b"noodles")?;

        let expected = [
            0xa0, 0x07, 0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x43, 0x04, 0x12, 0x05,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
