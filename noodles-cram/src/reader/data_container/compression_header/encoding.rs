use std::io;

use bytes::{Buf, Bytes};

use crate::{data_container::compression_header::Encoding, reader::num::get_itf8};

pub fn get_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let raw_kind = get_itf8(src)?;

    match raw_kind {
        0 => Ok(Encoding::Null),
        1 => get_external_encoding(src),
        2 => get_golomb_encoding(src),
        3 => get_huffman_encoding(src),
        4 => get_byte_array_len_encoding(src),
        5 => get_byte_array_stop_encoding(src),
        6 => get_beta_encoding(src),
        7 => get_subexp_encoding(src),
        8 => get_golomb_rice_encoding(src),
        9 => get_gamma_encoding(src),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid encoding kind",
        )),
    }
}

fn get_args(src: &mut Bytes) -> io::Result<Bytes> {
    let len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(src.split_to(len))
}

fn get_external_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;
    let block_content_id = get_itf8(&mut args)?;
    Ok(Encoding::External(block_content_id))
}

fn get_golomb_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let m = get_itf8(&mut args)?;

    Ok(Encoding::Golomb(offset, m))
}

fn get_byte_array_len_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let len_encoding = get_encoding(&mut args)?;
    let value_encoding = get_encoding(&mut args)?;

    Ok(Encoding::ByteArrayLen(
        Box::new(len_encoding),
        Box::new(value_encoding),
    ))
}

fn get_byte_array_stop_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    if !args.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let stop_byte = args.get_u8();
    let block_content_id = get_itf8(&mut args)?;

    Ok(Encoding::ByteArrayStop(stop_byte, block_content_id))
}

fn get_huffman_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let alphabet_len = get_itf8(&mut args).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut alphabet = Vec::with_capacity(alphabet_len);

    for _ in 0..alphabet_len {
        let symbol = get_itf8(&mut args)?;
        alphabet.push(symbol);
    }

    let bit_lens_len = get_itf8(&mut args).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bit_lens = Vec::with_capacity(bit_lens_len);

    for _ in 0..bit_lens_len {
        let len = get_itf8(&mut args).and_then(|n| {
            u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        bit_lens.push(len);
    }

    Ok(Encoding::Huffman(alphabet, bit_lens))
}

fn get_beta_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let len = get_itf8(&mut args).and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok(Encoding::Beta(offset, len))
}

fn get_subexp_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let k = get_itf8(&mut args)?;

    Ok(Encoding::Subexp(offset, k))
}

fn get_golomb_rice_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let log2_m = get_itf8(&mut args)?;

    Ok(Encoding::GolombRice(offset, log2_m))
}

fn get_gamma_encoding(src: &mut Bytes) -> io::Result<Encoding> {
    let mut args = get_args(src)?;
    let offset = get_itf8(&mut args)?;
    Ok(Encoding::Gamma(offset))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_null_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            0, // null encoding ID
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Null);

        Ok(())
    }

    #[test]
    fn test_get_external_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            1, // external encoding ID
            1, // args.len
            5, // block content ID
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::External(5));

        Ok(())
    }

    #[test]
    fn test_get_golomb_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2,  // Golomb encoding ID
            2,  // args.len
            1,  // offset
            10, // M
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Golomb(1, 10));

        Ok(())
    }

    #[test]
    fn test_get_huffman_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            3,  // Huffman encoding ID
            4,  // args.len
            1,  // alphabet.len
            65, // 'A'
            1,  // bit_lens.len
            0,  // 0
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Huffman(vec![65], vec![0]));

        Ok(())
    }

    #[test]
    fn test_get_byte_array_len_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            4,  // byte array len encoding ID
            6,  // args.len
            1,  // external encoding ID
            1,  // args.len
            13, // block content ID
            1,  // external encoding ID
            1,  // args.len
            21, // block content ID
        ]);

        let encoding = get_encoding(&mut data)?;

        assert_eq!(
            encoding,
            Encoding::ByteArrayLen(
                Box::new(Encoding::External(13)),
                Box::new(Encoding::External(21))
            )
        );

        Ok(())
    }

    #[test]
    fn test_get_byte_array_stop_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            5, // byte array stop encoding ID
            2, // args.len
            0, // NUL
            8, // block content ID
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::ByteArrayStop(0x00, 8));

        Ok(())
    }

    #[test]
    fn test_get_beta_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            6, // Beta encoding ID
            2, // args.len
            0, // offset
            8, // len
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Beta(0, 8));

        Ok(())
    }

    #[test]
    fn test_get_subexp_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            7, // subexponential encoding ID
            2, // args.len
            0, // offset
            1, // k
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Subexp(0, 1));

        Ok(())
    }

    #[test]
    fn test_get_golomb_rice_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            8, // Golomb-Rice encoding ID
            2, // args.len
            1, // offset
            3, // log2(M)
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::GolombRice(1, 3));

        Ok(())
    }

    #[test]
    fn test_get_gamma_encoding() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            9, // Elias gamma encoding ID
            1, // args.len
            1, // offset
        ]);

        let encoding = get_encoding(&mut data)?;
        assert_eq!(encoding, Encoding::Gamma(1));

        Ok(())
    }
}
