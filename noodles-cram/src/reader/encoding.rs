use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::{container::compression_header::Encoding, num::read_itf8};

pub fn read_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let raw_kind = read_itf8(reader)?;

    match raw_kind {
        0 => Ok(Encoding::Null),
        1 => read_external_encoding(reader),
        2 => unimplemented!("GOLOMB"),
        3 => read_huffman_encoding(reader),
        4 => read_byte_array_len_encoding(reader),
        5 => read_byte_array_stop_encoding(reader),
        6 => read_beta_encoding(reader),
        7 => unimplemented!("SUBEXP"),
        8 => unimplemented!("GOLOMB_RICE"),
        9 => unimplemented!("GAMMA"),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid encoding kind",
        )),
    }
}

fn read_args<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let len = read_itf8(reader)?;
    let mut buf = vec![0; len as usize];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn read_external_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let block_content_id = read_itf8(&mut args_reader)?;

    Ok(Encoding::External(block_content_id))
}

fn read_byte_array_len_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let len_encoding = read_encoding(&mut args_reader)?;
    let value_encoding = read_encoding(&mut args_reader)?;

    Ok(Encoding::ByteArrayLen(
        Box::new(len_encoding),
        Box::new(value_encoding),
    ))
}

fn read_byte_array_stop_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let stop_byte = args_reader.read_u8()?;
    let block_content_id = read_itf8(&mut args_reader)?;

    Ok(Encoding::ByteArrayStop(stop_byte, block_content_id))
}

fn read_huffman_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let alphabet_len = read_itf8(&mut args_reader)? as usize;
    let mut alphabet = Vec::with_capacity(alphabet_len);

    for _ in 0..alphabet_len {
        let symbol = read_itf8(&mut args_reader)?;
        alphabet.push(symbol);
    }

    let bit_lens_len = read_itf8(&mut args_reader)? as usize;
    let mut bit_lens = Vec::with_capacity(bit_lens_len);

    for _ in 0..bit_lens_len {
        let len = read_itf8(&mut args_reader)?;
        bit_lens.push(len);
    }

    Ok(Encoding::Huffman(alphabet, bit_lens))
}

fn read_beta_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader)?;
    let len = read_itf8(&mut args_reader)?;

    Ok(Encoding::Beta(offset, len))
}
