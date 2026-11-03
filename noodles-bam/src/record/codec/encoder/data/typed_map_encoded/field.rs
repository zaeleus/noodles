use std::{io, mem};

use memchr::memchr;
use noodles_sam::alignment::record::data::field::Tag;

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

const NUL: u8 = 0x00;

pub(super) fn write_field(dst: &mut Vec<u8>, src: &mut &[u8]) -> io::Result<()> {
    let original_src = *src;

    let tag = split_off_first_chunk::<2>(src).ok_or_else(unexpected_eof)?;
    let ty = src.split_off_first().copied().ok_or_else(unexpected_eof)?;

    if tag == Tag::CIGAR.as_ref() {
        return consume_cigar_field(src, ty);
    }

    match ty {
        CHARACTER_TYPE_VALUE | INT8_TYPE_VALUE | UINT8_TYPE_VALUE => {
            src.split_off_first().ok_or_else(unexpected_eof)?;
        }
        INT16_TYPE_VALUE | UINT16_TYPE_VALUE => {
            src.split_off(..mem::size_of::<u16>())
                .ok_or_else(unexpected_eof)?;
        }
        INT32_TYPE_VALUE | UINT32_TYPE_VALUE | FLOAT_TYPE_VALUE => {
            src.split_off(..mem::size_of::<u32>())
                .ok_or_else(unexpected_eof)?;
        }
        STRING_TYPE_VALUE => {
            let i = memchr(NUL, src).ok_or_else(unexpected_eof)?;

            // SAFETY: `i < src.len()`.
            let buf = src.split_off(..=i).unwrap();

            // SAFETY: `buf` is nonempty.
            let value = &buf[..buf.len() - 1];

            if !value.iter().all(|b| matches!(b, b' '..=b'~')) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid string",
                ));
            }
        }
        HEX_TYPE_VALUE => {
            let i = memchr(NUL, src).ok_or_else(unexpected_eof)?;

            // SAFETY: `i < src.len()`.
            let buf = src.split_off(..=i).unwrap();

            // SAFETY: `buf` is nonempty.
            let value = &buf[..buf.len() - 1];

            if !value.len().is_multiple_of(2)
                || !value.iter().all(|b| matches!(b, b'0'..=b'9' | b'A'..=b'F'))
            {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid hex"));
            }
        }
        ARRAY_TYPE_VALUE => {
            let subtype = src.split_off_first().copied().ok_or_else(unexpected_eof)?;

            let count = split_off_u32_le(src)
                .ok_or_else(unexpected_eof)
                .and_then(|n| {
                    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
                })?;

            let size = match subtype {
                INT8_TYPE_VALUE | UINT8_TYPE_VALUE => mem::size_of::<u8>(),
                INT16_TYPE_VALUE | UINT16_TYPE_VALUE => mem::size_of::<u16>(),
                INT32_TYPE_VALUE | UINT32_TYPE_VALUE | FLOAT_TYPE_VALUE => mem::size_of::<u32>(),
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "invalid subtype",
                    ));
                }
            };

            let len = size * count;
            src.split_off(..len).ok_or_else(unexpected_eof)?;
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid type")),
    }

    let len = original_src.len() - src.len();
    dst.extend(&original_src[..len]);

    Ok(())
}

fn consume_cigar_field(src: &mut &[u8], ty: u8) -> io::Result<()> {
    if ty != ARRAY_TYPE_VALUE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid CG field: expected B,I",
        ));
    }

    let subtype = src.split_off_first().copied().ok_or_else(unexpected_eof)?;

    if subtype != UINT32_TYPE_VALUE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid CG field: expected B,I",
        ));
    }

    let count = split_off_u32_le(src)
        .ok_or_else(unexpected_eof)
        .and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    let len = count * size_of::<u32>();
    src.split_off(..len).ok_or_else(unexpected_eof)?;

    Ok(())
}

fn split_off_first_chunk<'a, const N: usize>(src: &mut &'a [u8]) -> Option<&'a [u8; N]> {
    let (chunk, rest) = src.split_first_chunk()?;
    *src = rest;
    Some(chunk)
}

fn split_off_u32_le(src: &mut &[u8]) -> Option<u32> {
    split_off_first_chunk(src).map(|buf| u32::from_le_bytes(*buf))
}

fn unexpected_eof() -> io::Error {
    io::Error::from(io::ErrorKind::UnexpectedEof)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consume_cigar_field() -> io::Result<()> {
        let mut src = &[b'I', 0x01, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00][..];
        consume_cigar_field(&mut src, b'B')?;
        assert!(src.is_empty());

        assert!(matches!(
            consume_cigar_field(&mut &[0x01][..], b'C'),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        assert!(matches!(
            consume_cigar_field(&mut &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00][..], b'B'),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
