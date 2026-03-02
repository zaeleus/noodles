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
    file_definition::Version,
    io::reader::{
        collections::read_array,
        num::{read_signed_int, read_sint7_64, read_unsigned_int, read_unsigned_int_as},
    },
};

pub fn read_byte_encoding(src: &mut &[u8], version: Version) -> io::Result<Encoding<Byte>> {
    match read_kind(src, version)? {
        Kind::External => {
            let block_content_id = read_external_codec(src, version)?;
            Ok(Encoding::new(Byte::External { block_content_id }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = read_huffman_codec(src, version)?;
            Ok(Encoding::new(Byte::huffman(alphabet, bit_lens)))
        }
        Kind::ConstByte => {
            let value = read_const_byte_codec(src, version)?;
            Ok(Encoding::new(Byte::Constant { value }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<byte>: {kind:?}"),
        )),
    }
}

pub fn read_integer_encoding(src: &mut &[u8], version: Version) -> io::Result<Encoding<Integer>> {
    match read_kind(src, version)? {
        Kind::External => {
            let block_content_id = read_external_codec(src, version)?;
            Ok(Encoding::new(Integer::External { block_content_id }))
        }
        Kind::Golomb => {
            let (offset, m) = read_golomb_codec(src, version)?;
            Ok(Encoding::new(Integer::Golomb { offset, m }))
        }
        Kind::Huffman => {
            let (alphabet, bit_lens) = read_huffman_codec(src, version)?;
            Ok(Encoding::new(Integer::huffman(alphabet, bit_lens)))
        }
        Kind::Beta => {
            let (offset, len) = read_beta_codec(src, version)?;
            Ok(Encoding::new(Integer::Beta { offset, len }))
        }
        Kind::Subexp => {
            let (offset, k) = read_subexp_codec(src, version)?;
            Ok(Encoding::new(Integer::Subexp { offset, k }))
        }
        Kind::GolombRice => {
            let (offset, log2_m) = read_golomb_rice_codec(src, version)?;
            Ok(Encoding::new(Integer::GolombRice { offset, log2_m }))
        }
        Kind::Gamma => {
            let offset = read_gamma_codec(src, version)?;
            Ok(Encoding::new(Integer::Gamma { offset }))
        }
        Kind::VarintUnsigned => {
            let (block_content_id, offset) = read_varint_codec(src, version)?;
            Ok(Encoding::new(Integer::VarintUnsigned {
                block_content_id,
                offset,
            }))
        }
        Kind::VarintSigned => {
            let (block_content_id, offset) = read_varint_codec(src, version)?;
            Ok(Encoding::new(Integer::VarintSigned {
                block_content_id,
                offset,
            }))
        }
        Kind::ConstInt => {
            let value = read_const_int_codec(src, version)?;
            Ok(Encoding::new(Integer::ConstInt { value }))
        }
        kind => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid codec for encoding<int>: {kind:?}"),
        )),
    }
}

pub fn read_byte_array_encoding(
    src: &mut &[u8],
    version: Version,
) -> io::Result<Encoding<ByteArray>> {
    match read_kind(src, version)? {
        Kind::ByteArrayLength => {
            let (len_encoding, value_encoding) = read_byte_array_length_codec(src, version)?;

            Ok(Encoding::new(ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            }))
        }
        Kind::ByteArrayStop => {
            let (stop_byte, block_content_id) = read_byte_array_stop_codec(src, version)?;

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

/// Reads and discards a single encoding (of any type) from `src`.
///
/// The args block is length-prefixed and self-contained — even compound codecs like
/// `ByteArrayLength` store nested encodings within the args — so reading kind + args
/// is sufficient to skip the entire encoding.
pub fn consume_any_encoding(src: &mut &[u8], version: Version) -> io::Result<()> {
    read_kind(src, version)?;
    read_array(src, version)?;
    Ok(())
}

fn read_kind(src: &mut &[u8], version: Version) -> io::Result<Kind> {
    let n = read_unsigned_int(src, version)?;

    let kind = match n {
        0 => Kind::Null,
        1 => Kind::External,
        2 => Kind::Golomb,
        3 => Kind::Huffman,
        4 => Kind::ByteArrayLength,
        5 => Kind::ByteArrayStop,
        6 => Kind::Beta,
        7 => Kind::Subexp,
        8 => Kind::GolombRice,
        9 => Kind::Gamma,
        41 => Kind::VarintUnsigned,
        42 => Kind::VarintSigned,
        43 => Kind::ConstByte,
        44 => Kind::ConstInt,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid encoding kind: {n}"),
            ));
        }
    };

    if matches!(n, 41..=44) && version < Version::V4_0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("encoding kind {n} requires CRAM 4.0 or later"),
        ));
    }

    Ok(kind)
}

fn read_external_codec(src: &mut &[u8], version: Version) -> io::Result<block::ContentId> {
    let mut args = read_array(src, version)?;
    let block_content_id = read_unsigned_int(&mut args, version)?;
    Ok(block_content_id)
}

fn read_varint_codec(src: &mut &[u8], version: Version) -> io::Result<(block::ContentId, i64)> {
    let mut args = read_array(src, version)?;
    let block_content_id = read_unsigned_int(&mut args, version)?;
    let offset = read_sint7_64(&mut args)?;
    Ok((block_content_id, offset))
}

fn read_golomb_codec(src: &mut &[u8], version: Version) -> io::Result<(i32, i32)> {
    let mut args = read_array(src, version)?;

    let offset = read_signed_int(&mut args, version)?;
    let m = read_signed_int(&mut args, version)?;

    Ok((offset, m))
}

/// Reads a ByteArrayLength compound codec, which contains two nested sub-encodings
/// (a length encoding and a value encoding) serialized within its args block.
fn read_byte_array_length_codec(
    src: &mut &[u8],
    version: Version,
) -> io::Result<(Encoding<Integer>, Encoding<Byte>)> {
    let mut args = read_array(src, version)?;

    let len_encoding = read_integer_encoding(&mut args, version)?;
    let value_encoding = read_byte_encoding(&mut args, version)?;

    Ok((len_encoding, value_encoding))
}

fn read_byte_array_stop_codec(
    src: &mut &[u8],
    version: Version,
) -> io::Result<(u8, block::ContentId)> {
    let args = read_array(src, version)?;

    let (stop_byte, mut args) = args
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let block_content_id = read_unsigned_int(&mut args, version)?;

    Ok((*stop_byte, block_content_id))
}

fn read_huffman_codec(src: &mut &[u8], version: Version) -> io::Result<(Vec<i32>, Vec<u32>)> {
    let mut args = read_array(src, version)?;

    let alphabet_size: usize = read_unsigned_int_as(&mut args, version)?;

    let alphabet = (0..alphabet_size)
        .map(|_| read_signed_int(&mut args, version))
        .collect::<io::Result<_>>()?;

    let bit_lens_size: usize = read_unsigned_int_as(&mut args, version)?;

    let bit_lens = (0..bit_lens_size)
        .map(|_| read_unsigned_int_as::<_, u32>(&mut args, version))
        .collect::<io::Result<_>>()?;

    Ok((alphabet, bit_lens))
}

fn read_beta_codec(src: &mut &[u8], version: Version) -> io::Result<(i32, u32)> {
    let mut args = read_array(src, version)?;

    let offset = read_signed_int(&mut args, version)?;
    let len = read_unsigned_int_as(&mut args, version)?;

    Ok((offset, len))
}

fn read_subexp_codec(src: &mut &[u8], version: Version) -> io::Result<(i32, i32)> {
    let mut args = read_array(src, version)?;

    let offset = read_signed_int(&mut args, version)?;
    let k = read_signed_int(&mut args, version)?;

    Ok((offset, k))
}

fn read_golomb_rice_codec(src: &mut &[u8], version: Version) -> io::Result<(i32, i32)> {
    let mut args = read_array(src, version)?;

    let offset = read_signed_int(&mut args, version)?;
    let log2_m = read_signed_int(&mut args, version)?;

    Ok((offset, log2_m))
}

fn read_gamma_codec(src: &mut &[u8], version: Version) -> io::Result<i32> {
    let mut args = read_array(src, version)?;
    let offset = read_signed_int(&mut args, version)?;
    Ok(offset)
}

fn read_const_byte_codec(src: &mut &[u8], version: Version) -> io::Result<u8> {
    let args = read_array(src, version)?;

    let (&value, _) = args
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    Ok(value)
}

fn read_const_int_codec(src: &mut &[u8], version: Version) -> io::Result<i32> {
    let mut args = read_array(src, version)?;
    let value = read_signed_int(&mut args, version)?;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_kind() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Kind, version: Version) -> io::Result<()> {
            let actual = read_kind(&mut src, version)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        // Legacy kinds (valid for all versions)
        t(&[0x00], Kind::Null, Version::V3_0)?;
        t(&[0x01], Kind::External, Version::V3_0)?;
        t(&[0x02], Kind::Golomb, Version::V3_0)?;
        t(&[0x03], Kind::Huffman, Version::V3_0)?;
        t(&[0x04], Kind::ByteArrayLength, Version::V3_0)?;
        t(&[0x05], Kind::ByteArrayStop, Version::V3_0)?;
        t(&[0x06], Kind::Beta, Version::V3_0)?;
        t(&[0x07], Kind::Subexp, Version::V3_0)?;
        t(&[0x08], Kind::GolombRice, Version::V3_0)?;
        t(&[0x09], Kind::Gamma, Version::V3_0)?;

        // CRAM 4.0+ kinds
        t(&[0x29], Kind::VarintUnsigned, Version::V4_0)?;
        t(&[0x2a], Kind::VarintSigned, Version::V4_0)?;
        t(&[0x2b], Kind::ConstByte, Version::V4_0)?;
        t(&[0x2c], Kind::ConstInt, Version::V4_0)?;

        // Invalid kind
        assert!(matches!(
            read_kind(&mut &[0x0a][..], Version::V3_0),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        // CRAM 4.0 kinds rejected for pre-4.0 versions
        assert!(matches!(
            read_kind(&mut &[0x29][..], Version::V3_0),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[test]
    fn test_read_const_byte_codec() -> io::Result<()> {
        let src = [
            1,    // args.len
            0xff, // value
        ];

        let value = read_const_byte_codec(&mut &src[..], Version::V3_0)?;
        assert_eq!(value, 0xff);

        Ok(())
    }

    #[test]
    fn test_read_const_int_codec() -> io::Result<()> {
        let src = [
            1,  // args.len
            42, // value
        ];

        let value = read_const_int_codec(&mut &src[..], Version::V3_0)?;
        assert_eq!(value, 42);

        Ok(())
    }

    #[test]
    fn test_read_external_codec() -> io::Result<()> {
        let src = [
            1, // args.len
            5, // block content ID
        ];

        let block_content_id = read_external_codec(&mut &src[..], Version::V3_0)?;
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

        let (offset, m) = read_golomb_codec(&mut &src[..], Version::V3_0)?;
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

        let (alphabet, bit_lens) = read_huffman_codec(&mut &src[..], Version::V3_0)?;
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

        let (len_encoding, value_encoding) =
            read_byte_array_length_codec(&mut &src[..], Version::V3_0)?;

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

        let (stop_byte, block_content_id) =
            read_byte_array_stop_codec(&mut &src[..], Version::V3_0)?;
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

        let (offset, len) = read_beta_codec(&mut &src[..], Version::V3_0)?;
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

        let (offset, k) = read_subexp_codec(&mut &src[..], Version::V3_0)?;
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

        let (offset, log2_m) = read_golomb_rice_codec(&mut &src[..], Version::V3_0)?;
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

        let offset = read_gamma_codec(&mut &src[..], Version::V3_0)?;
        assert_eq!(offset, 1);

        Ok(())
    }

    #[test]
    fn test_read_varint_codec() -> io::Result<()> {
        let src = [
            2, // args.len
            5, // block content ID (uint7)
            0, // offset (sint7_64: 0)
        ];

        let (block_content_id, offset) = read_varint_codec(&mut &src[..], Version::V4_0)?;
        assert_eq!(block_content_id, 5);
        assert_eq!(offset, 0);

        Ok(())
    }

    #[test]
    fn test_read_integer_encoding_varint_unsigned() -> io::Result<()> {
        // VLQ-encoded: kind=41, args_len=2, block_content_id=3, offset=0
        let src = [41, 2, 3, 0];
        let encoding = read_integer_encoding(&mut &src[..], Version::V4_0)?;
        assert_eq!(
            encoding,
            Encoding::new(Integer::VarintUnsigned {
                block_content_id: 3,
                offset: 0,
            })
        );
        Ok(())
    }

    #[test]
    fn test_read_integer_encoding_varint_signed() -> io::Result<()> {
        // VLQ-encoded: kind=42, args_len=2, block_content_id=5, offset=0
        let src = [42, 2, 5, 0];
        let encoding = read_integer_encoding(&mut &src[..], Version::V4_0)?;
        assert_eq!(
            encoding,
            Encoding::new(Integer::VarintSigned {
                block_content_id: 5,
                offset: 0,
            })
        );
        Ok(())
    }

    #[test]
    fn test_read_integer_encoding_const_int() -> io::Result<()> {
        // VLQ-encoded: kind=44, args_len=1, value=84 (zigzag-encoded 42)
        let src = [44, 1, 84];
        let encoding = read_integer_encoding(&mut &src[..], Version::V4_0)?;
        assert_eq!(encoding, Encoding::new(Integer::ConstInt { value: 42 }));
        Ok(())
    }

    #[test]
    fn test_read_byte_encoding_const_byte_v4() -> io::Result<()> {
        // VLQ-encoded: kind=43, args_len=1, value=0x41
        let src = [43, 1, 0x41];
        let encoding = read_byte_encoding(&mut &src[..], Version::V4_0)?;
        assert_eq!(encoding, Encoding::new(Byte::Constant { value: 0x41 }));
        Ok(())
    }
}
