use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::{data_container::compression_header::Encoding, reader::num::read_itf8};

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
        7 => read_subexp_encoding(reader),
        8 => unimplemented!("GOLOMB_RICE"),
        9 => read_gamma_encoding(reader),
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

fn read_subexp_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader)?;
    let k = read_itf8(&mut args_reader)?;

    Ok(Encoding::Subexp(offset, k))
}

fn read_gamma_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let args = read_args(reader)?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader)?;

    Ok(Encoding::Gamma(offset))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_null_encoding() -> io::Result<()> {
        let data = [
            0, // null encoding ID
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::Null);

        Ok(())
    }

    #[test]
    fn test_read_external_encoding() -> io::Result<()> {
        let data = [
            1, // external encoding ID
            1, // args.len
            5, // block content ID
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::External(5));

        Ok(())
    }

    #[test]
    fn test_read_huffman_encoding() -> io::Result<()> {
        let data = [
            3,  // Huffman encoding ID
            4,  // args.len
            1,  // alphabet.len
            65, // 'A'
            1,  // bit_lens.len
            0,  // 0
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::Huffman(vec![65], vec![0]));

        Ok(())
    }

    #[test]
    fn test_read_byte_array_len_encoding() -> io::Result<()> {
        let data = [
            4,  // byte array len encoding ID
            6,  // args.len
            1,  // external encoding ID
            1,  // args.len
            13, // block content ID
            1,  // external encoding ID
            1,  // args.len
            21, // block content ID
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;

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
    fn test_read_byte_array_stop_encoding() -> io::Result<()> {
        let data = [
            5, // byte array stop encoding ID
            2, // args.len
            0, // NUL
            8, // block content ID
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::ByteArrayStop(0x00, 8));

        Ok(())
    }

    #[test]
    fn test_read_beta_encoding() -> io::Result<()> {
        let data = [
            6, // Beta encoding ID
            2, // args.len
            0, // offset
            8, // len
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::Beta(0, 8));

        Ok(())
    }

    #[test]
    fn test_read_subexp_encoding() -> io::Result<()> {
        let data = [
            7, // subexponential encoding ID
            2, // args.len
            0, // offset
            1, // k
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::Subexp(0, 1));

        Ok(())
    }

    #[test]
    fn test_read_gamma_encoding() -> io::Result<()> {
        let data = [
            9, // Elias gamma encoding ID
            1, // args.len
            1, // offset
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader)?;
        assert_eq!(encoding, Encoding::Gamma(1));

        Ok(())
    }
}
