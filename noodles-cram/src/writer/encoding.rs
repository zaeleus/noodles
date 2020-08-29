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
        Encoding::Huffman(..) => todo!("HUFFMAN"),
        Encoding::ByteArrayLen(len_encoding, value_encoding) => {
            write_byte_array_len_encoding(writer, len_encoding, value_encoding)
        }
        Encoding::ByteArrayStop(stop_byte, block_content_id) => {
            write_byte_array_stop_encoding(writer, *stop_byte, *block_content_id)
        }
        Encoding::Beta(..) => todo!("BETA"),
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
