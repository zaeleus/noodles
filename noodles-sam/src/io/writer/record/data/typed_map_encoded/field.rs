use std::io::{self, Write};

use bstr::ByteSlice;

use crate::io::writer::num;

const FIELD_SEPARATOR: u8 = b'\t';
const COMPONENT_SEPARATOR: u8 = b':';
const ARRAY_ITEM_SEPARATOR: u8 = b',';

const CHARACTER_TYPE_VALUE: u8 = b'A';
const INT8_TYPE_VALUE: u8 = b'c';
const UINT8_TYPE_VALUE: u8 = b'C';
const INT16_TYPE_VALUE: u8 = b's';
const UINT16_TYPE_VALUE: u8 = b'S';
const INT32_TYPE_VALUE: u8 = b'i';
const UINT32_TYPE_VALUE: u8 = b'I';
const FLOAT_TYPE_VALUE: u8 = b'f';
const STRING_TYPE_VALUE: u8 = b'Z';
const HEX_TYPE_VALUE: u8 = b'H';
const ARRAY_TYPE_VALUE: u8 = b'B';

const INTEGER_TYPE_VALUE: u8 = INT32_TYPE_VALUE;

const NUL: u8 = 0x00;

type Tag = [u8; 2];

pub(super) fn write_field<W>(writer: &mut W, src: &mut &[u8]) -> io::Result<()>
where
    W: Write,
{
    let tag = read_tag(src)?;
    let ty = read_u8(src)?;

    match ty {
        CHARACTER_TYPE_VALUE => write_character_field(writer, src, tag)?,
        INT8_TYPE_VALUE => write_i8_field(writer, src, tag)?,
        UINT8_TYPE_VALUE => write_u8_field(writer, src, tag)?,
        INT16_TYPE_VALUE => write_i16_field(writer, src, tag)?,
        UINT16_TYPE_VALUE => write_u16_field(writer, src, tag)?,
        INT32_TYPE_VALUE => write_i32_field(writer, src, tag)?,
        UINT32_TYPE_VALUE => write_u32_field(writer, src, tag)?,
        FLOAT_TYPE_VALUE => write_f32_field(writer, src, tag)?,
        STRING_TYPE_VALUE => write_string_field(writer, src, tag)?,
        HEX_TYPE_VALUE => write_hex_field(writer, src, tag)?,
        ARRAY_TYPE_VALUE => write_array_field(writer, src, tag)?,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid type")),
    }

    Ok(())
}

fn read_tag(src: &mut &[u8]) -> io::Result<Tag> {
    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;

    // § 1.5 "The alignment section: optional fields" (2025-08-12): "`[A-Za-z][A-Za-z0-9]`".
    if buf[0].is_ascii_alphabetic() && buf[1].is_ascii_alphanumeric() {
        Ok(*buf)
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
}

fn write_character_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, CHARACTER_TYPE_VALUE)?;

    let n = read_u8(src)?;

    if !n.is_ascii_graphic() {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    writer.write_all(&[n])?;

    Ok(())
}

fn write_i8_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let n = read_u8(src)?;
    num::write_i8(writer, n as i8)?;

    Ok(())
}

fn write_u8_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let n = read_u8(src)?;
    num::write_u8(writer, n)?;

    Ok(())
}

fn write_i16_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;
    let n = i16::from_le_bytes(*buf);
    num::write_i16(writer, n)?;

    Ok(())
}

fn write_u16_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;
    let n = u16::from_le_bytes(*buf);
    num::write_u16(writer, n)?;

    Ok(())
}

fn write_i32_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;
    let n = i32::from_le_bytes(*buf);
    num::write_i32(writer, n)?;

    Ok(())
}

fn write_u32_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, INTEGER_TYPE_VALUE)?;

    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;
    let n = u32::from_le_bytes(*buf);
    num::write_u32(writer, n)?;

    Ok(())
}

fn write_f32_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, FLOAT_TYPE_VALUE)?;

    let buf = split_off_first_chunk(src).ok_or_else(unexpected_eof)?;
    let n = f32::from_le_bytes(*buf);

    if !n.is_finite() {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    num::write_f32(writer, n)?;

    Ok(())
}

fn write_string_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, STRING_TYPE_VALUE)?;

    let i = src.as_bstr().find_byte(NUL).ok_or_else(unexpected_eof)?;

    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];

    if !buf.iter().all(|b| matches!(b, b' '..=b'~')) {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    writer.write_all(buf)?;

    Ok(())
}

fn write_hex_field<W>(writer: &mut W, src: &mut &[u8], tag: Tag) -> io::Result<()>
where
    W: Write,
{
    write_field_prefix(writer, tag, HEX_TYPE_VALUE)?;

    let i = src.as_bstr().find_byte(NUL).ok_or_else(unexpected_eof)?;

    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];

    if !buf.len().is_multiple_of(2) || !buf.iter().all(|b| matches!(b, b'0'..=b'9' | b'A'..=b'F')) {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    writer.write_all(buf)?;

    Ok(())
}

fn write_field_prefix<W>(writer: &mut W, tag: Tag, ty: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[
        FIELD_SEPARATOR,
        tag[0],
        tag[1],
        COMPONENT_SEPARATOR,
        ty,
        COMPONENT_SEPARATOR,
    ])
}

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

    let len = count;
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

    let len = count;
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
    split_off_first_chunk(src)
        .map(|buf| u32::from_le_bytes(*buf))
        .ok_or_else(unexpected_eof)
}

fn split_off_first_chunk<'a, const N: usize>(src: &mut &'a [u8]) -> Option<&'a [u8; N]> {
    let (chunk, rest) = src.split_first_chunk()?;
    *src = rest;
    Some(chunk)
}

fn unexpected_eof() -> io::Error {
    io::Error::from(io::ErrorKind::UnexpectedEof)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut src: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_field(buf, &mut src)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ZAAn", b"\tZA:A:n")?;
        t(&mut buf, &[b'Z', b'B', b'c', 0x00], b"\tZB:i:0")?;
        t(&mut buf, &[b'Z', b'C', b'C', 0x00], b"\tZC:i:0")?;
        t(&mut buf, &[b'Z', b'D', b's', 0x00, 0x00], b"\tZD:i:0")?;
        t(&mut buf, &[b'Z', b'E', b'S', 0x00, 0x00], b"\tZE:i:0")?;
        t(
            &mut buf,
            &[b'Z', b'F', b'i', 0x00, 0x00, 0x00, 0x00],
            b"\tZF:i:0",
        )?;
        t(
            &mut buf,
            &[b'Z', b'G', b'I', 0x00, 0x00, 0x00, 0x00],
            b"\tZG:i:0",
        )?;
        t(
            &mut buf,
            &[b'Z', b'H', b'f', 0x00, 0x00, 0x00, 0x00],
            b"\tZH:f:0",
        )?;
        t(
            &mut buf,
            &[b'Z', b'I', b'Z', b'n', b'd', b'l', b's', 0x00],
            b"\tZI:Z:ndls",
        )?;
        t(
            &mut buf,
            &[b'Z', b'J', b'H', b'C', b'A', b'F', b'E', 0x00],
            b"\tZJ:H:CAFE",
        )?;
        t(
            &mut buf,
            &[b'Z', b'K', b'B', b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            b"\tZK:B:c,0",
        )?;

        Ok(())
    }

    #[test]
    fn test_read_tag() -> io::Result<()> {
        assert_eq!(read_tag(&mut &b"NH"[..])?, *b"NH");

        assert!(matches!(
            read_tag(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_tag(&mut &b"n\t"[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_character_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let tag = *b"ZZ";

        buf.clear();
        write_character_field(&mut buf, &mut &b"n"[..], tag)?;
        assert_eq!(buf, b"\tZZ:A:n");

        buf.clear();
        assert!(matches!(
            write_character_field(&mut buf, &mut &[][..], tag),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        assert!(matches!(
            write_character_field(&mut buf, &mut &b"\n"[..], tag),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_f32_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let tag = *b"ZZ";

        buf.clear();
        write_f32_field(&mut buf, &mut &[0x00, 0x00, 0x00, 0x00][..], tag)?;
        assert_eq!(buf, b"\tZZ:f:0");

        buf.clear();
        assert!(matches!(
            write_f32_field(&mut buf, &mut &[][..], tag),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        let src = f32::NAN.to_le_bytes();
        assert!(matches!(
            write_f32_field(&mut buf, &mut &src[..], tag),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_string_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let tag = *b"ZZ";

        buf.clear();
        let src = [b'n', b'd', b'l', b's', 0x00];
        write_string_field(&mut buf, &mut &src[..], tag)?;
        assert_eq!(buf, b"\tZZ:Z:ndls");

        buf.clear();
        assert!(matches!(
            write_string_field(&mut buf, &mut &[][..], tag),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        let src = [b'n', b'd', b'\t', b'l', b's', 0x00];
        assert!(matches!(
            write_string_field(&mut buf, &mut &src[..], tag),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_hex_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let tag = *b"ZZ";

        buf.clear();
        let src = [b'C', b'A', b'F', b'E', 0x00];
        write_hex_field(&mut buf, &mut &src[..], tag)?;
        assert_eq!(buf, b"\tZZ:H:CAFE");

        buf.clear();
        assert!(matches!(
            write_hex_field(&mut buf, &mut &[][..], tag),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        let src = [b'n', b'd', b'l', b's', 0x00];
        assert!(matches!(
            write_hex_field(&mut buf, &mut &src[..], tag),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

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
