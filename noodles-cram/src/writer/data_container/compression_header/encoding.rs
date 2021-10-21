use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::{
    data_container::compression_header::{encoding, Encoding},
    num::Itf8,
    writer::num::write_itf8,
};

pub fn write_encoding<W>(writer: &mut W, encoding: &Encoding) -> io::Result<()>
where
    W: Write,
{
    match encoding {
        Encoding::Null => write_null_encoding(writer),
        Encoding::External(block_content_id) => write_external_encoding(writer, *block_content_id),
        Encoding::Golomb(..) => unimplemented!("GOLOMB"),
        Encoding::Huffman(alphabet, bit_lens) => write_huffman_encoding(writer, alphabet, bit_lens),
        Encoding::ByteArrayLen(len_encoding, value_encoding) => {
            write_byte_array_len_encoding(writer, len_encoding, value_encoding)
        }
        Encoding::ByteArrayStop(stop_byte, block_content_id) => {
            write_byte_array_stop_encoding(writer, *stop_byte, *block_content_id)
        }
        Encoding::Beta(offset, len) => write_beta_encoding(writer, *offset, *len),
        Encoding::Subexp(offset, k) => write_subexp_encoding(writer, *offset, *k),
        Encoding::GolombRice(..) => unimplemented!("GOLOMB_RICE"),
        Encoding::Gamma(offset) => write_gamma_encoding(writer, *offset),
    }
}

fn write_args<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len = buf.len() as Itf8;
    write_itf8(writer, len)?;
    writer.write_all(buf)
}

fn write_null_encoding<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    write_itf8(writer, i32::from(encoding::Kind::Null))
}

fn write_external_encoding<W>(writer: &mut W, block_content_id: Itf8) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, block_content_id)?;

    write_itf8(writer, i32::from(encoding::Kind::External))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_huffman_encoding<W>(writer: &mut W, alphabet: &[i32], bit_lens: &[u32]) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    let alphabet_len = alphabet.len() as Itf8;
    write_itf8(&mut args, alphabet_len)?;

    for &symbol in alphabet {
        write_itf8(&mut args, symbol)?;
    }

    let bit_lens_len = bit_lens.len() as Itf8;
    write_itf8(&mut args, bit_lens_len)?;

    for &len in bit_lens {
        let len =
            Itf8::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_itf8(&mut args, len)?;
    }

    write_itf8(writer, i32::from(encoding::Kind::Huffman))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_byte_array_len_encoding<W>(
    writer: &mut W,
    len_encoding: &Encoding,
    value_encoding: &Encoding,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();

    write_encoding(&mut args, len_encoding)?;
    write_encoding(&mut args, value_encoding)?;

    write_itf8(writer, i32::from(encoding::Kind::ByteArrayLen))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_byte_array_stop_encoding<W>(
    writer: &mut W,
    stop_byte: u8,
    block_content_id: Itf8,
) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    args.write_u8(stop_byte)?;
    write_itf8(&mut args, block_content_id)?;

    write_itf8(writer, i32::from(encoding::Kind::ByteArrayStop))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_beta_encoding<W>(writer: &mut W, offset: Itf8, len: u32) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;

    let len = Itf8::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut args, len)?;

    write_itf8(writer, i32::from(encoding::Kind::Beta))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_subexp_encoding<W>(writer: &mut W, offset: Itf8, k: Itf8) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;
    write_itf8(&mut args, k)?;

    write_itf8(writer, i32::from(encoding::Kind::Subexp))?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_gamma_encoding<W>(writer: &mut W, offset: Itf8) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;

    write_itf8(writer, i32::from(encoding::Kind::Gamma))?;
    write_args(writer, &args)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_null_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_null_encoding(&mut buf)?;

        let expected = [
            0, // null encoding ID
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_external_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_external_encoding(&mut buf, 5)?;

        let expected = [
            1, // external encoding ID
            1, // args.len
            5, // block content ID
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_huffman_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_huffman_encoding(&mut buf, &[65], &[0])?;

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
    fn test_write_byte_array_len_encoding() -> io::Result<()> {
        let mut buf = Vec::new();

        let len_encoding = Encoding::External(13);
        let value_encoding = Encoding::External(21);

        write_byte_array_len_encoding(&mut buf, &len_encoding, &value_encoding)?;

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
    fn test_write_byte_array_stop_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_byte_array_stop_encoding(&mut buf, 0x00, 8)?;

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
    fn test_write_beta_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_beta_encoding(&mut buf, 0, 8)?;

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
    fn test_write_subexp_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_subexp_encoding(&mut buf, 0, 1)?;

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
    fn test_write_gamma_encoding() -> io::Result<()> {
        let mut buf = Vec::new();
        write_gamma_encoding(&mut buf, 1)?;

        let expected = [
            9, // Elias gamma encoding ID
            1, // args.len
            1, // offset
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
