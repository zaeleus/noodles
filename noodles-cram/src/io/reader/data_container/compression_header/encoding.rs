use std::io;

use bytes::{Buf, Bytes};

use crate::{
    container::block,
    data_container::compression_header::{
        encoding::{
            codec::{Byte, ByteArray, Integer},
            Kind,
        },
        Encoding,
    },
    io::reader::num::get_itf8,
};

pub fn get_encoding_for_byte_codec(src: &mut Bytes) -> io::Result<Encoding<Byte>> {
    match get_kind(src)? {
        Kind::External => {
            let block_content_id = get_external_codec(src)?;
            Ok(Encoding::new(Byte::External { block_content_id }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = get_huffman_codec(src)?;
            Ok(Encoding::new(Byte::Huffman { alphabet, bit_lens }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<byte>: {kind:?}"),
        )),
    }
}

pub fn get_encoding_for_integer_codec(src: &mut Bytes) -> io::Result<Encoding<Integer>> {
    match get_kind(src)? {
        Kind::External => {
            let block_content_id = get_external_codec(src)?;
            Ok(Encoding::new(Integer::External { block_content_id }))
        }
        Kind::Golomb => {
            let (offset, m) = get_golomb_codec(src)?;
            Ok(Encoding::new(Integer::Golomb { offset, m }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = get_huffman_codec(src)?;
            Ok(Encoding::new(Integer::Huffman { alphabet, bit_lens }))
        }
        Kind::Beta => {
            let (offset, len) = get_beta_codec(src)?;
            Ok(Encoding::new(Integer::Beta { offset, len }))
        }
        Kind::Subexp => {
            let (offset, k) = get_subexp_codec(src)?;
            Ok(Encoding::new(Integer::Subexp { offset, k }))
        }
        Kind::GolombRice => {
            let (offset, log2_m) = get_golomb_rice_codec(src)?;
            Ok(Encoding::new(Integer::GolombRice { offset, log2_m }))
        }
        Kind::Gamma => {
            let offset = get_gamma_codec(src)?;
            Ok(Encoding::new(Integer::Gamma { offset }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<int>: {kind:?}"),
        )),
    }
}

pub fn get_encoding_for_byte_array_codec(src: &mut Bytes) -> io::Result<Encoding<ByteArray>> {
    match get_kind(src)? {
        Kind::ByteArrayLen => {
            let (len_encoding, value_encoding) = get_byte_array_len_codec(src)?;

            Ok(Encoding::new(ByteArray::ByteArrayLen {
                len_encoding,
                value_encoding,
            }))
        }
        Kind::ByteArrayStop => {
            let (stop_byte, block_content_id) = get_byte_array_stop_codec(src)?;

            Ok(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte,
                block_content_id,
            }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<byte[]>: {kind:?}"),
        )),
    }
}

fn get_kind(src: &mut Bytes) -> io::Result<Kind> {
    match get_itf8(src)? {
        0 => Ok(Kind::Null),
        1 => Ok(Kind::External),
        2 => Ok(Kind::Golomb),
        3 => Ok(Kind::Huffman),
        4 => Ok(Kind::ByteArrayLen),
        5 => Ok(Kind::ByteArrayStop),
        6 => Ok(Kind::Beta),
        7 => Ok(Kind::Subexp),
        8 => Ok(Kind::GolombRice),
        9 => Ok(Kind::Gamma),
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

fn get_external_codec(src: &mut Bytes) -> io::Result<block::ContentId> {
    let mut args = get_args(src)?;
    let block_content_id = get_itf8(&mut args).map(block::ContentId::from)?;
    Ok(block_content_id)
}

fn get_golomb_codec(src: &mut Bytes) -> io::Result<(i32, i32)> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let m = get_itf8(&mut args)?;

    Ok((offset, m))
}

fn get_byte_array_len_codec(src: &mut Bytes) -> io::Result<(Encoding<Integer>, Encoding<Byte>)> {
    let mut args = get_args(src)?;

    let len_encoding = get_encoding_for_integer_codec(&mut args)?;
    let value_encoding = get_encoding_for_byte_codec(&mut args)?;

    Ok((len_encoding, value_encoding))
}

fn get_byte_array_stop_codec(src: &mut Bytes) -> io::Result<(u8, block::ContentId)> {
    let mut args = get_args(src)?;

    let stop_byte = args
        .try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))?;

    let block_content_id = get_itf8(&mut args).map(block::ContentId::from)?;

    Ok((stop_byte, block_content_id))
}

fn get_huffman_codec(src: &mut Bytes) -> io::Result<(Vec<i32>, Vec<u32>)> {
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

    Ok((alphabet, bit_lens))
}

fn get_beta_codec(src: &mut Bytes) -> io::Result<(i32, u32)> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let len = get_itf8(&mut args).and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok((offset, len))
}

fn get_subexp_codec(src: &mut Bytes) -> io::Result<(i32, i32)> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let k = get_itf8(&mut args)?;

    Ok((offset, k))
}

fn get_golomb_rice_codec(src: &mut Bytes) -> io::Result<(i32, i32)> {
    let mut args = get_args(src)?;

    let offset = get_itf8(&mut args)?;
    let log2_m = get_itf8(&mut args)?;

    Ok((offset, log2_m))
}

fn get_gamma_codec(src: &mut Bytes) -> io::Result<i32> {
    let mut args = get_args(src)?;
    let offset = get_itf8(&mut args)?;
    Ok(offset)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_kind() -> io::Result<()> {
        fn t(buf: &'static [u8], expected: Kind) -> io::Result<()> {
            let mut src = Bytes::from_static(buf);
            let actual = get_kind(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], Kind::Null)?;
        t(&[0x01], Kind::External)?;
        t(&[0x02], Kind::Golomb)?;
        t(&[0x03], Kind::Huffman)?;
        t(&[0x04], Kind::ByteArrayLen)?;
        t(&[0x05], Kind::ByteArrayStop)?;
        t(&[0x06], Kind::Beta)?;
        t(&[0x07], Kind::Subexp)?;
        t(&[0x08], Kind::GolombRice)?;
        t(&[0x09], Kind::Gamma)?;

        let mut src = Bytes::from_static(&[0x0a]);
        assert!(matches!(
            get_kind(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[test]
    fn test_get_external_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            1, // args.len
            5, // block content ID
        ]);

        let block_content_id = get_external_codec(&mut data)?;
        assert_eq!(block_content_id, block::ContentId::from(5));

        Ok(())
    }

    #[test]
    fn test_get_golomb_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2,  // args.len
            1,  // offset
            10, // M
        ]);

        let (offset, m) = get_golomb_codec(&mut data)?;
        assert_eq!(offset, 1);
        assert_eq!(m, 10);

        Ok(())
    }

    #[test]
    fn test_get_huffman_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            4,  // args.len
            1,  // alphabet.len
            65, // 'A'
            1,  // bit_lens.len
            0,  // 0
        ]);

        let (alphabet, bit_lens) = get_huffman_codec(&mut data)?;
        assert_eq!(alphabet, [65]);
        assert_eq!(bit_lens, [0]);

        Ok(())
    }

    #[test]
    fn test_get_byte_array_len_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            6,  // args.len
            1,  // external encoding ID
            1,  // args.len
            13, // block content ID
            1,  // external encoding ID
            1,  // args.len
            21, // block content ID
        ]);

        let (len_encoding, value_encoding) = get_byte_array_len_codec(&mut data)?;

        assert_eq!(
            len_encoding,
            Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(13)
            })
        );
        assert_eq!(
            value_encoding,
            Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(21)
            })
        );

        Ok(())
    }

    #[test]
    fn test_get_byte_array_stop_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2, // args.len
            0, // NUL
            8, // block content ID
        ]);

        let (stop_byte, block_content_id) = get_byte_array_stop_codec(&mut data)?;
        assert_eq!(stop_byte, 0);
        assert_eq!(block_content_id, block::ContentId::from(8));

        Ok(())
    }

    #[test]
    fn test_get_beta_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2, // args.len
            0, // offset
            8, // len
        ]);

        let (offset, len) = get_beta_codec(&mut data)?;
        assert_eq!(offset, 0);
        assert_eq!(len, 8);

        Ok(())
    }

    #[test]
    fn test_get_subexp_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2, // args.len
            0, // offset
            1, // k
        ]);

        let (offset, k) = get_subexp_codec(&mut data)?;
        assert_eq!(offset, 0);
        assert_eq!(k, 1);

        Ok(())
    }

    #[test]
    fn test_get_golomb_rice_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            2, // args.len
            1, // offset
            3, // log2(M)
        ]);

        let (offset, log2_m) = get_golomb_rice_codec(&mut data)?;
        assert_eq!(offset, 1);
        assert_eq!(log2_m, 3);

        Ok(())
    }

    #[test]
    fn test_get_gamma_codec() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            1, // args.len
            1, // offset
        ]);

        let offset = get_gamma_codec(&mut data)?;
        assert_eq!(offset, 1);

        Ok(())
    }
}
