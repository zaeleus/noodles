use std::io::{self, Write};

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
    io::writer::num::{write_int, write_signed_int, write_u8},
};

pub fn write_byte_encoding<W>(
    writer: &mut W,
    encoding: &Encoding<Byte>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        Byte::External { block_content_id } => {
            write_external_codec(writer, *block_content_id, version)
        }
        Byte::Huffman {
            alphabet, bit_lens, ..
        } => write_huffman_codec(writer, alphabet, bit_lens, version),
        Byte::Constant { value } => write_const_byte_codec(writer, *value, version),
    }
}

pub fn write_integer_encoding<W>(
    writer: &mut W,
    encoding: &Encoding<Integer>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        Integer::External { block_content_id } => {
            write_external_codec(writer, *block_content_id, version)
        }
        Integer::Golomb { offset, m } => write_golomb_codec(writer, *offset, *m, version),
        Integer::Huffman {
            alphabet, bit_lens, ..
        } => write_huffman_codec(writer, alphabet, bit_lens, version),
        Integer::Beta { offset, len } => write_beta_codec(writer, *offset, *len, version),
        Integer::Subexp { offset, k } => write_subexp_codec(writer, *offset, *k, version),
        Integer::GolombRice { offset, log2_m } => {
            write_golomb_rice_codec(writer, *offset, *log2_m, version)
        }
        Integer::Gamma { offset } => write_gamma_codec(writer, *offset, version),
        Integer::VarintUnsigned {
            block_content_id,
            offset,
        } => write_varint_unsigned_codec(writer, *block_content_id, *offset, version),
        Integer::VarintSigned {
            block_content_id,
            offset,
        } => write_varint_signed_codec(writer, *block_content_id, *offset, version),
        Integer::ConstInt { value } => write_const_int_codec(writer, *value, version),
    }
}

pub fn write_byte_array_encoding<W>(
    writer: &mut W,
    encoding: &Encoding<ByteArray>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        ByteArray::ByteArrayLength {
            len_encoding,
            value_encoding,
        } => write_byte_array_length_codec(writer, len_encoding, value_encoding, version),
        ByteArray::ByteArrayStop {
            stop_byte,
            block_content_id,
        } => write_byte_array_stop_codec(writer, *stop_byte, *block_content_id, version),
    }
}

fn write_kind<W>(writer: &mut W, kind: Kind, version: Version) -> io::Result<()>
where
    W: Write,
{
    // Legacy kinds 0–9 are valid for all versions; kinds 41–44 are CRAM 4.0+ only.
    debug_assert!(
        !matches!(
            kind,
            Kind::VarintUnsigned | Kind::VarintSigned | Kind::ConstByte | Kind::ConstInt
        ) || version >= Version::V4_0,
        "CRAM 4.0+ codec kind {kind:?} used with version {version:?}"
    );

    let n = match kind {
        Kind::Null => 0,
        Kind::External => 1,
        Kind::Golomb => 2,
        Kind::Huffman => 3,
        Kind::ByteArrayLength => 4,
        Kind::ByteArrayStop => 5,
        Kind::Beta => 6,
        Kind::Subexp => 7,
        Kind::GolombRice => 8,
        Kind::Gamma => 9,
        Kind::VarintUnsigned => 41,
        Kind::VarintSigned => 42,
        Kind::ConstByte => 43,
        Kind::ConstInt => 44,
    };

    write_int(writer, version, n)
}

fn write_args<W>(writer: &mut W, buf: &[u8], version: Version) -> io::Result<()>
where
    W: Write,
{
    let len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, len)?;
    writer.write_all(buf)
}

/// Writes unsigned codec args (block content IDs, array lengths, bit lengths).
fn write_unsigned_arg<W>(buf: &mut W, version: Version, value: i32) -> io::Result<()>
where
    W: Write,
{
    write_int(buf, version, value)
}

/// Writes signed codec args (offsets, alphabet symbols, codec parameters).
///
/// Uses zigzag encoding for CRAM 4.0, matching the reader's `read_signed_int`.
fn write_signed_arg<W>(buf: &mut W, version: Version, value: i32) -> io::Result<()>
where
    W: Write,
{
    write_signed_int(buf, version, value)
}

fn write_external_codec<W>(
    writer: &mut W,
    block_content_id: block::ContentId,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    write_unsigned_arg(&mut args, version, block_content_id)?;

    write_kind(writer, Kind::External, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_golomb_codec<W>(writer: &mut W, offset: i32, m: i32, version: Version) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, offset)?;
    write_signed_arg(&mut args, version, m)?;

    write_kind(writer, Kind::Golomb, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_huffman_codec<W>(
    writer: &mut W,
    alphabet: &[i32],
    bit_lens: &[u32],
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    let alphabet_len = i32::try_from(alphabet.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_unsigned_arg(&mut args, version, alphabet_len)?;

    for &symbol in alphabet {
        write_signed_arg(&mut args, version, symbol)?;
    }

    let bit_lens_len = i32::try_from(bit_lens.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_unsigned_arg(&mut args, version, bit_lens_len)?;

    for &len in bit_lens {
        let len = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_unsigned_arg(&mut args, version, len)?;
    }

    write_kind(writer, Kind::Huffman, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_byte_array_length_codec<W>(
    writer: &mut W,
    len_encoding: &Encoding<Integer>,
    value_encoding: &Encoding<Byte>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    write_integer_encoding(&mut args, len_encoding, version)?;
    write_byte_encoding(&mut args, value_encoding, version)?;

    write_kind(writer, Kind::ByteArrayLength, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_byte_array_stop_codec<W>(
    writer: &mut W,
    stop_byte: u8,
    block_content_id: block::ContentId,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_u8(&mut args, stop_byte)?;

    write_unsigned_arg(&mut args, version, block_content_id)?;

    write_kind(writer, Kind::ByteArrayStop, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_beta_codec<W>(writer: &mut W, offset: i32, len: u32, version: Version) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, offset)?;

    let len = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_unsigned_arg(&mut args, version, len)?;

    write_kind(writer, Kind::Beta, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_subexp_codec<W>(writer: &mut W, offset: i32, k: i32, version: Version) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, offset)?;
    write_signed_arg(&mut args, version, k)?;

    write_kind(writer, Kind::Subexp, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_golomb_rice_codec<W>(
    writer: &mut W,
    offset: i32,
    log2_m: i32,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, offset)?;
    write_signed_arg(&mut args, version, log2_m)?;

    write_kind(writer, Kind::GolombRice, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_gamma_codec<W>(writer: &mut W, offset: i32, version: Version) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, offset)?;

    write_kind(writer, Kind::Gamma, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_varint_unsigned_codec<W>(
    writer: &mut W,
    block_content_id: block::ContentId,
    offset: i64,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::num::write_sint7_64;

    let mut args = Vec::new();
    write_unsigned_arg(&mut args, version, block_content_id)?;
    write_sint7_64(&mut args, offset)?;

    write_kind(writer, Kind::VarintUnsigned, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_varint_signed_codec<W>(
    writer: &mut W,
    block_content_id: block::ContentId,
    offset: i64,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::num::write_sint7_64;

    let mut args = Vec::new();
    write_unsigned_arg(&mut args, version, block_content_id)?;
    write_sint7_64(&mut args, offset)?;

    write_kind(writer, Kind::VarintSigned, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_const_byte_codec<W>(writer: &mut W, value: u8, version: Version) -> io::Result<()>
where
    W: Write,
{
    let args = [value];

    write_kind(writer, Kind::ConstByte, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

fn write_const_int_codec<W>(writer: &mut W, value: i32, version: Version) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_signed_arg(&mut args, version, value)?;

    write_kind(writer, Kind::ConstInt, version)?;
    write_args(writer, &args, version)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_external_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_external_codec(&mut buf, 5, Version::default())?;

        let expected = [
            1, // external encoding ID
            1, // args.len
            5, // block content ID
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_golomb_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_golomb_codec(&mut buf, 1, 10, Version::default())?;

        let expected = [
            2,  // Golomb encoding ID
            2,  // args.len
            1,  // offset
            10, // m
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_huffman_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_huffman_codec(&mut buf, &[65], &[0], Version::default())?;

        let expected = [
            3,  // Huffman encoding ID
            4,  // args.len
            1,  // alphabet.len
            65, // 'A'
            1,  // bit_lens.len
            0,  // 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_byte_array_length_codec() -> io::Result<()> {
        let mut buf = Vec::new();

        let len_encoding = Encoding::new(Integer::External {
            block_content_id: 13,
        });
        let value_encoding = Encoding::new(Byte::External {
            block_content_id: 21,
        });

        write_byte_array_length_codec(
            &mut buf,
            &len_encoding,
            &value_encoding,
            Version::default(),
        )?;

        let expected = [
            4,  // byte array len encoding ID
            6,  // args.len
            1,  // external encoding ID
            1,  // args.len
            13, // block content ID
            1,  // external encoding ID
            1,  // args.len
            21, // block content ID
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_byte_array_stop_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_byte_array_stop_codec(&mut buf, 0x00, 8, Version::default())?;

        let expected = [
            5, // byte array stop encoding ID
            2, // args.len
            0, // NUL
            8, // block content ID
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_beta_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_beta_codec(&mut buf, 0, 8, Version::default())?;

        let expected = [
            6, // Beta encoding ID
            2, // args.len
            0, // offset
            8, // len
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_subexp_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_subexp_codec(&mut buf, 0, 1, Version::default())?;

        let expected = [
            7, // subexponential encoding ID
            2, // args.len
            0, // offset
            1, // k
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_golomb_rice_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_golomb_rice_codec(&mut buf, 1, 3, Version::default())?;

        let expected = [
            8, // Golomb encoding ID
            2, // args.len
            1, // offset
            3, // m
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_gamma_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_gamma_codec(&mut buf, 1, Version::default())?;

        let expected = [
            9, // Elias gamma encoding ID
            1, // args.len
            1, // offset
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_varint_unsigned_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_varint_unsigned_codec(&mut buf, 5, 0, Version::V4_0)?;

        let expected = [
            41, // varint unsigned encoding ID
            2,  // args.len
            5,  // block content ID
            0,  // offset (sint7_64: 0)
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_varint_signed_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_varint_signed_codec(&mut buf, 3, 0, Version::V4_0)?;

        let expected = [
            42, // varint signed encoding ID
            2,  // args.len
            3,  // block content ID
            0,  // offset (sint7_64: 0)
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_const_byte_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_const_byte_codec(&mut buf, 0x41, Version::V4_0)?;

        let expected = [
            43,   // const byte encoding ID
            1,    // args.len
            0x41, // value = 'A'
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_const_int_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_const_int_codec(&mut buf, 42, Version::V4_0)?;

        let expected = [
            44, // const int encoding ID
            1,  // args.len
            84, // value (zigzag-encoded: 42 * 2 = 84)
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
