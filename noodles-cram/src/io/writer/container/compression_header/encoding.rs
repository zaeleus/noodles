use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::{
    container::{
        block,
        compression_header::{
            encoding::{
                codec::{Byte, ByteArray, Integer},
                Kind,
            },
            Encoding,
        },
    },
    io::writer::num::write_itf8,
};

pub fn write_encoding_for_byte_codec<W>(writer: &mut W, encoding: &Encoding<Byte>) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        Byte::External { block_content_id } => write_external_codec(writer, *block_content_id),
        Byte::Huffman { alphabet, bit_lens } => write_huffman_codec(writer, alphabet, bit_lens),
    }
}

pub fn write_encoding_for_integer_codec<W>(
    writer: &mut W,
    encoding: &Encoding<Integer>,
) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        Integer::External { block_content_id } => write_external_codec(writer, *block_content_id),
        Integer::Golomb { offset, m } => write_golomb_codec(writer, *offset, *m),
        Integer::Huffman { alphabet, bit_lens } => write_huffman_codec(writer, alphabet, bit_lens),
        Integer::Beta { offset, len } => write_beta_codec(writer, *offset, *len),
        Integer::Subexp { offset, k } => write_subexp_codec(writer, *offset, *k),
        Integer::GolombRice { offset, log2_m } => write_golomb_rice_codec(writer, *offset, *log2_m),
        Integer::Gamma { offset } => write_gamma_codec(writer, *offset),
    }
}

pub fn write_encoding_for_byte_array_codec<W>(
    writer: &mut W,
    encoding: &Encoding<ByteArray>,
) -> io::Result<()>
where
    W: Write,
{
    match encoding.get() {
        ByteArray::ByteArrayLen {
            len_encoding,
            value_encoding,
        } => write_byte_array_len_codec(writer, len_encoding, value_encoding),
        ByteArray::ByteArrayStop {
            stop_byte,
            block_content_id,
        } => write_byte_array_stop_codec(writer, *stop_byte, *block_content_id),
    }
}

fn write_kind<W>(writer: &mut W, kind: Kind) -> io::Result<()>
where
    W: Write,
{
    let n = match kind {
        Kind::Null => 0,
        Kind::External => 1,
        Kind::Golomb => 2,
        Kind::Huffman => 3,
        Kind::ByteArrayLen => 4,
        Kind::ByteArrayStop => 5,
        Kind::Beta => 6,
        Kind::Subexp => 7,
        Kind::GolombRice => 8,
        Kind::Gamma => 9,
    };

    write_itf8(writer, n)
}

fn write_args<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;
    writer.write_all(buf)
}

fn write_external_codec<W>(writer: &mut W, block_content_id: block::ContentId) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    write_itf8(&mut args, block_content_id)?;

    write_kind(writer, Kind::External)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_golomb_codec<W>(writer: &mut W, offset: i32, m: i32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;
    write_itf8(&mut args, m)?;

    write_kind(writer, Kind::Golomb)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_huffman_codec<W>(writer: &mut W, alphabet: &[i32], bit_lens: &[u32]) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    let alphabet_len = i32::try_from(alphabet.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut args, alphabet_len)?;

    for &symbol in alphabet {
        write_itf8(&mut args, symbol)?;
    }

    let bit_lens_len = i32::try_from(bit_lens.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut args, bit_lens_len)?;

    for &len in bit_lens {
        let len = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_itf8(&mut args, len)?;
    }

    write_kind(writer, Kind::Huffman)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_byte_array_len_codec<W>(
    writer: &mut W,
    len_encoding: &Encoding<Integer>,
    value_encoding: &Encoding<Byte>,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    write_encoding_for_integer_codec(&mut args, len_encoding)?;
    write_encoding_for_byte_codec(&mut args, value_encoding)?;

    write_kind(writer, Kind::ByteArrayLen)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_byte_array_stop_codec<W>(
    writer: &mut W,
    stop_byte: u8,
    block_content_id: block::ContentId,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    args.write_u8(stop_byte)?;

    write_itf8(&mut args, block_content_id)?;

    write_kind(writer, Kind::ByteArrayStop)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_beta_codec<W>(writer: &mut W, offset: i32, len: u32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;

    let len = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut args, len)?;

    write_kind(writer, Kind::Beta)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_subexp_codec<W>(writer: &mut W, offset: i32, k: i32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;
    write_itf8(&mut args, k)?;

    write_kind(writer, Kind::Subexp)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_golomb_rice_codec<W>(writer: &mut W, offset: i32, log2_m: i32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;
    write_itf8(&mut args, log2_m)?;

    write_kind(writer, Kind::GolombRice)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_gamma_codec<W>(writer: &mut W, offset: i32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;

    write_kind(writer, Kind::Gamma)?;
    write_args(writer, &args)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_external_codec() -> io::Result<()> {
        let mut buf = Vec::new();
        write_external_codec(&mut buf, 5)?;

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
        write_golomb_codec(&mut buf, 1, 10)?;

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
        write_huffman_codec(&mut buf, &[65], &[0])?;

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
    fn test_write_byte_array_len_codec() -> io::Result<()> {
        let mut buf = Vec::new();

        let len_encoding = Encoding::new(Integer::External {
            block_content_id: 13,
        });
        let value_encoding = Encoding::new(Byte::External {
            block_content_id: 21,
        });

        write_byte_array_len_codec(&mut buf, &len_encoding, &value_encoding)?;

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
        write_byte_array_stop_codec(&mut buf, 0x00, 8)?;

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
        write_beta_codec(&mut buf, 0, 8)?;

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
        write_subexp_codec(&mut buf, 0, 1)?;

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
        write_golomb_rice_codec(&mut buf, 1, 3)?;

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
        write_gamma_codec(&mut buf, 1)?;

        let expected = [
            9, // Elias gamma encoding ID
            1, // args.len
            1, // offset
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
