//! Writes BAM-encoded data fields as SAM text.
//!
//! This transcodes the binary layout directly. Going through
//! [`crate::alignment::record::data::field::Value`] instead means materializing a value for every
//! field of every record, and reaching it through a boxed iterator.

use std::io::{self, Write};

use bstr::BStr;

use super::field::value::{write_float, write_hex, write_string};
use crate::io::writer::num;

const DELIMITER: u8 = b':';
const SEPARATOR: u8 = b',';
const NUL: u8 = 0x00;

pub(super) fn write_field_encoded_data<W>(writer: &mut W, mut src: &[u8]) -> io::Result<()>
where
    W: Write,
{
    const TAB: u8 = b'\t';

    while !src.is_empty() {
        writer.write_all(&[TAB])?;
        write_field(writer, &mut src)?;
    }

    Ok(())
}

fn write_field<W>(writer: &mut W, src: &mut &[u8]) -> io::Result<()>
where
    W: Write,
{
    let tag = split_off::<2>(src)?;
    writer.write_all(&tag)?;
    writer.write_all(&[DELIMITER])?;

    let ty = split_off_first(src)?;

    // § 1.5 "The alignment section: optional fields" (2024-11-06): the integer types collapse to
    // `i` in SAM.
    let sam_ty = match ty {
        b'A' => b'A',
        b'c' | b'C' | b's' | b'S' | b'i' | b'I' => b'i',
        b'f' => b'f',
        b'Z' => b'Z',
        b'H' => b'H',
        b'B' => b'B',
        _ => return Err(invalid_data("invalid data field type")),
    };

    writer.write_all(&[sam_ty])?;
    writer.write_all(&[DELIMITER])?;

    write_value(writer, src, ty)
}

fn write_value<W>(writer: &mut W, src: &mut &[u8], ty: u8) -> io::Result<()>
where
    W: Write,
{
    match ty {
        b'A' => writer.write_all(&[split_off_first(src)?]),
        b'c' => num::write_i8(writer, split_off_first(src)? as i8),
        b'C' => num::write_u8(writer, split_off_first(src)?),
        b's' => num::write_i16(writer, i16::from_le_bytes(split_off(src)?)),
        b'S' => num::write_u16(writer, u16::from_le_bytes(split_off(src)?)),
        b'i' => num::write_i32(writer, i32::from_le_bytes(split_off(src)?)),
        b'I' => num::write_u32(writer, u32::from_le_bytes(split_off(src)?)),
        b'f' => write_float(writer, f32::from_le_bytes(split_off(src)?)),
        b'Z' => write_string(writer, split_off_nul_terminated(src)?),
        b'H' => write_hex(writer, split_off_nul_terminated(src)?),
        b'B' => write_array(writer, src),
        _ => Err(invalid_data("invalid data field type")),
    }
}

fn write_array<W>(writer: &mut W, src: &mut &[u8]) -> io::Result<()>
where
    W: Write,
{
    let subtype = split_off_first(src)?;

    if !matches!(subtype, b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f') {
        return Err(invalid_data("invalid array subtype"));
    }

    writer.write_all(&[subtype])?;

    let len = u32::from_le_bytes(split_off(src)?);

    for _ in 0..len {
        writer.write_all(&[SEPARATOR])?;
        write_value(writer, src, subtype)?;
    }

    Ok(())
}

fn split_off<const N: usize>(src: &mut &[u8]) -> io::Result<[u8; N]> {
    if src.len() < N {
        return Err(unexpected_eof());
    }

    let (buf, rest) = src.split_at(N);
    *src = rest;

    // SAFETY: `buf` is `N` bytes.
    Ok(buf.try_into().unwrap())
}

fn split_off_first(src: &mut &[u8]) -> io::Result<u8> {
    split_off::<1>(src).map(|buf| buf[0])
}

fn split_off_nul_terminated<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    let i = src
        .iter()
        .position(|&b| b == NUL)
        .ok_or_else(unexpected_eof)?;
    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];
    Ok(BStr::new(buf))
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn unexpected_eof() -> io::Error {
    io::Error::from(io::ErrorKind::UnexpectedEof)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field_encoded_data() -> io::Result<()> {
        fn t(src: &[u8], expected: &[u8]) -> io::Result<()> {
            let mut buf = Vec::new();
            write_field_encoded_data(&mut buf, src)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        t(&[], b"")?;
        t(&[b'N', b'H', b'C', 0x01], b"\tNH:i:1")?;
        t(&[b'N', b'H', b'c', 0xff], b"\tNH:i:-1")?;
        t(&[b'N', b'H', b's', 0xff, 0xff], b"\tNH:i:-1")?;
        t(&[b'N', b'H', b'S', 0x01, 0x00], b"\tNH:i:1")?;
        t(&[b'N', b'H', b'i', 0xff, 0xff, 0xff, 0xff], b"\tNH:i:-1")?;
        t(&[b'N', b'H', b'I', 0x01, 0x00, 0x00, 0x00], b"\tNH:i:1")?;
        t(b"COAn", b"\tCO:A:n")?;
        t(&[b'F', b'Z', b'f', 0x00, 0x00, 0x80, 0x3f], b"\tFZ:f:1")?;
        t(
            &[b'C', b'O', b'Z', b'n', b'd', b'l', b's', 0x00],
            b"\tCO:Z:ndls",
        )?;
        t(
            &[b'C', b'O', b'H', b'C', b'A', b'F', b'E', 0x00],
            b"\tCO:H:CAFE",
        )?;

        t(
            &[b'F', b'Z', b'B', b'C', 0x02, 0x00, 0x00, 0x00, 0x07, 0x0d],
            b"\tFZ:B:C,7,13",
        )?;

        // Two fields in one buffer.
        t(
            &[
                b'N', b'H', b'C', 0x01, b'R', b'G', b'Z', b'r', b'g', b'0', 0x00,
            ],
            b"\tNH:i:1\tRG:Z:rg0",
        )?;

        Ok(())
    }

    #[test]
    fn test_write_field_encoded_data_with_invalid_input() {
        // A truncated value.
        let mut buf = Vec::new();
        assert!(matches!(
            write_field_encoded_data(&mut buf, &[b'N', b'H', b'i', 0x00]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        // An unknown type.
        let mut buf = Vec::new();
        assert!(matches!(
            write_field_encoded_data(&mut buf, &[b'N', b'H', b'?', 0x00]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        // An unterminated string.
        let mut buf = Vec::new();
        assert!(matches!(
            write_field_encoded_data(&mut buf, b"COZn"),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
    }
}
