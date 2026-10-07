#![allow(dead_code)]

use std::io::{self, Write};

use crate::io::writer::num;

const FIELD_SEPARATOR: u8 = b'\t';
const COMPONENT_SEPARATOR: u8 = b':';
const ARRAY_ITEM_SEPARATOR: u8 = b',';

const INT8_TYPE_VALUE: u8 = b'c';
const UINT8_TYPE_VALUE: u8 = b'C';
const INT16_TYPE_VALUE: u8 = b's';
const UINT16_TYPE_VALUE: u8 = b'S';
const INT32_TYPE_VALUE: u8 = b'i';
const UINT32_TYPE_VALUE: u8 = b'I';
const FLOAT_TYPE_VALUE: u8 = b'f';
const ARRAY_TYPE_VALUE: u8 = b'B';

type Tag = [u8; 2];

fn write_array_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    let subtype = read_u8(src)?;

    let count = read_u32_le(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    match subtype {
        INT8_TYPE_VALUE => write_i8_array_field(writer, src, tag, count)?,
        UINT8_TYPE_VALUE => write_u8_array_field(writer, src, tag, count)?,
        INT16_TYPE_VALUE => write_i16_array_field(writer, src, tag, count)?,
        UINT16_TYPE_VALUE => write_u16_array_field(writer, src, tag, count)?,
        INT32_TYPE_VALUE => write_i32_array_field(writer, src, tag, count)?,
        UINT32_TYPE_VALUE => write_u32_array_field(writer, src, tag, count)?,
        FLOAT_TYPE_VALUE => write_f32_array_field(writer, src, tag, count)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid subtype",
            ));
        }
    }

    Ok(())
}

fn write_i8_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, INT8_TYPE_VALUE)?;

    let len = count * size_of::<i8>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    for &n in buf {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;
        num::write_i8(writer, n as i8)?;
    }

    Ok(())
}

fn write_u8_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, UINT8_TYPE_VALUE)?;

    let len = count * size_of::<u8>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    for &n in buf {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;
        num::write_u8(writer, n)?;
    }

    Ok(())
}

fn write_i16_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, INT16_TYPE_VALUE)?;

    let len = count * size_of::<i16>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    let (chunks, []) = buf.as_chunks() else {
        unreachable!();
    };

    for chunk in chunks {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;

        let n = i16::from_le_bytes(*chunk);
        num::write_i16(writer, n)?;
    }

    Ok(())
}

fn write_u16_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, UINT16_TYPE_VALUE)?;

    let len = count * size_of::<u16>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    let (chunks, []) = buf.as_chunks() else {
        unreachable!();
    };

    for chunk in chunks {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;

        let n = u16::from_le_bytes(*chunk);
        num::write_u16(writer, n)?;
    }

    Ok(())
}

fn write_i32_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, INT32_TYPE_VALUE)?;

    let len = count * size_of::<i32>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    let (chunks, []) = buf.as_chunks() else {
        unreachable!();
    };

    for chunk in chunks {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;

        let n = i32::from_le_bytes(*chunk);
        num::write_i32(writer, n)?;
    }

    Ok(())
}

fn write_u32_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, UINT32_TYPE_VALUE)?;

    let len = count * size_of::<u32>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    let (chunks, []) = buf.as_chunks() else {
        unreachable!();
    };

    for chunk in chunks {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;

        let n = u32::from_le_bytes(*chunk);
        num::write_u32(writer, n)?;
    }

    Ok(())
}

fn write_f32_array_field<W>(
    writer: &mut W,
    src: &mut &[u8],
    tag: Tag,
    count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_array_field_prefix(writer, tag, FLOAT_TYPE_VALUE)?;

    let len = count * size_of::<f32>();
    let buf = src.split_off(..len).ok_or_else(unexpected_eof)?;

    let (chunks, []) = buf.as_chunks() else {
        unreachable!();
    };

    for chunk in chunks {
        writer.write_all(&[ARRAY_ITEM_SEPARATOR])?;

        let n = f32::from_le_bytes(*chunk);

        if !n.is_finite() {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        num::write_f32(writer, n)?;
    }

    Ok(())
}

fn write_array_field_prefix<W>(writer: &mut W, tag: Tag, subtype: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[
        FIELD_SEPARATOR,
        tag[0],
        tag[1],
        COMPONENT_SEPARATOR,
        ARRAY_TYPE_VALUE,
        COMPONENT_SEPARATOR,
        subtype,
    ])
}

fn read_u8(src: &mut &[u8]) -> io::Result<u8> {
    src.split_off_first().copied().ok_or_else(unexpected_eof)
}

fn read_u32_le(src: &mut &[u8]) -> io::Result<u32> {
    let (buf, rest) = src.split_first_chunk().ok_or_else(unexpected_eof)?;
    *src = rest;
    Ok(u32::from_le_bytes(*buf))
}

fn unexpected_eof() -> io::Error {
    io::Error::from(io::ErrorKind::UnexpectedEof)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_array_field() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut src: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_array_field(buf, &mut src, *b"ZZ")?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:c,0",
        )?;
        t(
            &mut buf,
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:C,0",
        )?;
        t(
            &mut buf,
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:s,0",
        )?;
        t(
            &mut buf,
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:S,0",
        )?;
        t(
            &mut buf,
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:i,0",
        )?;
        t(
            &mut buf,
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:I,0",
        )?;
        t(
            &mut buf,
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            b"\tZZ:B:f,0",
        )?;

        Ok(())
    }

    #[test]
    fn test_write_f32_array_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let tag = *b"ZZ";

        buf.clear();
        let src = &[0x00, 0x00, 0x00, 0x00];
        write_f32_array_field(&mut buf, &mut &src[..], tag, 1)?;
        assert_eq!(buf, b"\tZZ:B:f,0");

        buf.clear();
        assert!(matches!(
            write_f32_array_field(&mut buf, &mut &[][..], tag, 1),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        let src = f32::NAN.to_le_bytes();
        assert!(matches!(
            write_f32_array_field(&mut buf, &mut &src[..], tag, 1),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
