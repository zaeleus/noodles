use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::{
    container::compression_header::Encoding,
    num::{write_itf8, Itf8},
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
        Encoding::Subexp(..) => unimplemented!("SUBEXP"),
        Encoding::GolombRice(..) => unimplemented!("GOLOMB_RICE"),
        Encoding::Gamma(_) => unimplemented!("GAMMA"),
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
    // TODO: convert from encoding
    write_itf8(writer, 0)
}

fn write_external_encoding<W>(writer: &mut W, block_content_id: Itf8) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, block_content_id)?;

    // TODO: convert from encoding
    write_itf8(writer, 1)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_huffman_encoding<W>(writer: &mut W, alphabet: &[i32], bit_lens: &[i32]) -> io::Result<()>
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
        write_itf8(&mut args, len)?;
    }

    // TODO: convert from encoding
    write_itf8(writer, 3)?;
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

    // TODO: convert from encoding
    write_itf8(writer, 4)?;
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

    // TODO: convert from encoding
    write_itf8(writer, 5)?;
    write_args(writer, &args)?;

    Ok(())
}

fn write_beta_encoding<W>(writer: &mut W, offset: Itf8, len: Itf8) -> io::Result<()>
where
    W: Write,
{
    let mut args = Vec::new();
    write_itf8(&mut args, offset)?;
    write_itf8(&mut args, len)?;

    // TODO: convert from encoding
    write_itf8(writer, 6)?;
    write_args(writer, &args)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
