use std::io;

use crate::{
    container::{
        block,
        compression_header::{
            Encoding,
            encoding::{
                Kind,
                codec::{Byte, ByteArray, Integer},
            },
        },
    },
    io::reader::{
        collections::read_array,
        num::{read_itf8, read_itf8_as},
    },
};

pub fn read_byte_encoding(src: &mut &[u8]) -> io::Result<Encoding<Byte>> {
    match read_kind(src)? {
        Kind::External => {
            let block_content_id = read_external_codec(src)?;
            Ok(Encoding::new(Byte::External { block_content_id }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = read_huffman_codec(src)?;
            Ok(Encoding::new(Byte::Huffman { alphabet, bit_lens }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<byte>: {kind:?}"),
        )),
    }
}

pub fn read_integer_encoding(src: &mut &[u8]) -> io::Result<Encoding<Integer>> {
    match read_kind(src)? {
        Kind::External => {
            let block_content_id = read_external_codec(src)?;
            Ok(Encoding::new(Integer::External { block_content_id }))
        }
        Kind::Golomb => {
            let (offset, m) = read_golomb_codec(src)?;
            Ok(Encoding::new(Integer::Golomb { offset, m }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = read_huffman_codec(src)?;
            Ok(Encoding::new(Integer::Huffman { alphabet, bit_lens }))
        }
        Kind::Beta => {
            let (offset, len) = read_beta_codec(src)?;
            Ok(Encoding::new(Integer::Beta { offset, len }))
        }
        Kind::Subexp => {
            let (offset, k) = read_subexp_codec(src)?;
            Ok(Encoding::new(Integer::Subexp { offset, k }))
        }
        Kind::GolombRice => {
            let (offset, log2_m) = read_golomb_rice_codec(src)?;
            Ok(Encoding::new(Integer::GolombRice { offset, log2_m }))
        }
        Kind::Gamma => {
            let offset = read_gamma_codec(src)?;
            Ok(Encoding::new(Integer::Gamma { offset }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<int>: {kind:?}"),
        )),
    }
}

pub fn read_byte_array_encoding(src: &mut &[u8]) -> io::Result<Encoding<ByteArray>> {
    match read_kind(src)? {
        Kind::ByteArrayLength => {
            let (len_encoding, value_encoding) = read_byte_array_length_codec(src)?;

            Ok(Encoding::new(ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            }))
        }
        Kind::ByteArrayStop => {
            let (stop_byte, block_content_id) = read_byte_array_stop_codec(src)?;

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

pub fn consume_any_encoding(src: &mut &[u8]) -> io::Result<()> {
    read_kind(src)?;
    read_array(src)?;
    Ok(())
}

fn read_kind(src: &mut &[u8]) -> io::Result<Kind> {
    match read_itf8(src)? {
        0 => Ok(Kind::Null),
        1 => Ok(Kind::External),
        2 => Ok(Kind::Golomb),
        3 => Ok(Kind::Huffman),
        4 => Ok(Kind::ByteArrayLength),
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

fn read_external_codec(src: &mut &[u8]) -> io::Result<block::ContentId> {
    let mut args = read_array(src)?;
    let block_content_id = read_itf8(&mut args)?;
    Ok(block_content_id)
}

fn read_golomb_codec(src: &mut &[u8]) -> io::Result<(i32, i32)> {
    let mut args = read_array(src)?;

    let offset = read_itf8(&mut args)?;
    let m = read_itf8(&mut args)?;

    Ok((offset, m))
}

fn read_byte_array_length_codec(
    src: &mut &[u8],
) -> io::Result<(Encoding<Integer>, Encoding<Byte>)> {
    let mut args = read_array(src)?;

    let len_encoding = read_integer_encoding(&mut args)?;
    let value_encoding = read_byte_encoding(&mut args)?;

    Ok((len_encoding, value_encoding))
}

fn read_byte_array_stop_codec(src: &mut &[u8]) -> io::Result<(u8, block::ContentId)> {
    let args = read_array(src)?;

    let (stop_byte, mut args) = args
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let block_content_id = read_itf8(&mut args)?;

    Ok((*stop_byte, block_content_id))
}

fn read_huffman_codec(src: &mut &[u8]) -> io::Result<(Vec<i32>, Vec<u32>)> {
    let mut args = read_array(src)?;

    let alphabet_size: usize = read_itf8_as(&mut args)?;

    let alphabet = (0..alphabet_size)
        .map(|_| read_itf8(&mut args))
        .collect::<io::Result<_>>()?;

    let bit_lens_size: usize = read_itf8_as(&mut args)?;

    let bit_lens = (0..bit_lens_size)
        .map(|_| {
            read_itf8(&mut args).and_then(|n| {
                u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
        })
        .collect::<io::Result<_>>()?;

    Ok((alphabet, bit_lens))
}

fn read_beta_codec(src: &mut &[u8]) -> io::Result<(i32, u32)> {
    let mut args = read_array(src)?;

    let offset = read_itf8(&mut args)?;
    let len = read_itf8_as(&mut args)?;

    Ok((offset, len))
}

fn read_subexp_codec(src: &mut &[u8]) -> io::Result<(i32, i32)> {
    let mut args = read_array(src)?;

    let offset = read_itf8(&mut args)?;
    let k = read_itf8(&mut args)?;

    Ok((offset, k))
}

fn read_golomb_rice_codec(src: &mut &[u8]) -> io::Result<(i32, i32)> {
    let mut args = read_array(src)?;

    let offset = read_itf8(&mut args)?;
    let log2_m = read_itf8(&mut args)?;

    Ok((offset, log2_m))
}

fn read_gamma_codec(src: &mut &[u8]) -> io::Result<i32> {
    let mut args = read_array(src)?;
    let offset = read_itf8(&mut args)?;
    Ok(offset)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_kind() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Kind) -> io::Result<()> {
            let actual = read_kind(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], Kind::Null)?;
        t(&[0x01], Kind::External)?;
        t(&[0x02], Kind::Golomb)?;
        t(&[0x03], Kind::Huffman)?;
        t(&[0x04], Kind::ByteArrayLength)?;
        t(&[0x05], Kind::ByteArrayStop)?;
        t(&[0x06], Kind::Beta)?;
        t(&[0x07], Kind::Subexp)?;
        t(&[0x08], Kind::GolombRice)?;
        t(&[0x09], Kind::Gamma)?;

        assert!(matches!(
            read_kind(&mut &[0x0a][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[test]
    fn test_read_external_codec() -> io::Result<()> {
        let src = [
            1, // args.len
            5, // block content ID
        ];

        let block_content_id = read_external_codec(&mut &src[..])?;
        assert_eq!(block_content_id, 5);

        Ok(())
    }

    #[test]
    fn test_read_golomb_codec() -> io::Result<()> {
        let src = [
            2,  // args.len
            1,  // offset
            10, // M
        ];

        let (offset, m) = read_golomb_codec(&mut &src[..])?;
        assert_eq!(offset, 1);
        assert_eq!(m, 10);

        Ok(())
    }

    #[test]
    fn test_read_huffman_codec() -> io::Result<()> {
        let src = [
            4,  // args.len
            1,  // alphabet.len
            65, // 'A'
            1,  // bit_lens.len
            0,  // 0
        ];

        let (alphabet, bit_lens) = read_huffman_codec(&mut &src[..])?;
        assert_eq!(alphabet, [65]);
        assert_eq!(bit_lens, [0]);

        Ok(())
    }

    #[test]
    fn test_read_byte_array_length_codec() -> io::Result<()> {
        let src = [
            6,  // args.len
            1,  // external encoding ID
            1,  // args.len
            13, // block content ID
            1,  // external encoding ID
            1,  // args.len
            21, // block content ID
        ];

        let (len_encoding, value_encoding) = read_byte_array_length_codec(&mut &src[..])?;

        assert_eq!(
            len_encoding,
            Encoding::new(Integer::External {
                block_content_id: 13
            })
        );
        assert_eq!(
            value_encoding,
            Encoding::new(Byte::External {
                block_content_id: 21
            })
        );

        Ok(())
    }

    #[test]
    fn test_read_byte_array_stop_codec() -> io::Result<()> {
        let src = [
            2, // args.len
            0, // NUL
            8, // block content ID
        ];

        let (stop_byte, block_content_id) = read_byte_array_stop_codec(&mut &src[..])?;
        assert_eq!(stop_byte, 0);
        assert_eq!(block_content_id, 8);

        Ok(())
    }

    #[test]
    fn test_read_beta_codec() -> io::Result<()> {
        let src = [
            2, // args.len
            0, // offset
            8, // len
        ];

        let (offset, len) = read_beta_codec(&mut &src[..])?;
        assert_eq!(offset, 0);
        assert_eq!(len, 8);

        Ok(())
    }

    #[test]
    fn test_read_subexp_codec() -> io::Result<()> {
        let src = [
            2, // args.len
            0, // offset
            1, // k
        ];

        let (offset, k) = read_subexp_codec(&mut &src[..])?;
        assert_eq!(offset, 0);
        assert_eq!(k, 1);

        Ok(())
    }

    #[test]
    fn test_read_golomb_rice_codec() -> io::Result<()> {
        let src = [
            2, // args.len
            1, // offset
            3, // log2(M)
        ];

        let (offset, log2_m) = read_golomb_rice_codec(&mut &src[..])?;
        assert_eq!(offset, 1);
        assert_eq!(log2_m, 3);

        Ok(())
    }

    #[test]
    fn test_read_gamma_codec() -> io::Result<()> {
        let src = [
            1, // args.len
            1, // offset
        ];

        let offset = read_gamma_codec(&mut &src[..])?;
        assert_eq!(offset, 1);

        Ok(())
    }
}
