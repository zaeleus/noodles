//! BAM record data field writer.

pub mod field;

use std::{io, mem};

use memchr::memchr;
use noodles_sam::alignment::record::{Data, DataRef, data::field::Tag};

use self::field::write_field;

const NUL: u8 = 0x00;

pub(super) fn write_data(dst: &mut Vec<u8>, data: DataRef<'_>) -> io::Result<()> {
    match data {
        DataRef::FieldEncoded(src) => write_field_encoded_data(dst, src),
        DataRef::Data(d) => write_generic_data(dst, d),
    }
}

fn write_field_encoded_data(dst: &mut Vec<u8>, src: &[u8]) -> io::Result<()> {
    if is_valid(src)? {
        dst.extend(src);
        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
}

fn write_generic_data<D>(dst: &mut Vec<u8>, data: D) -> io::Result<()>
where
    D: Data,
{
    for result in data.iter() {
        let (tag, value) = result?;

        if &tag == Tag::CIGAR.as_ref() {
            continue;
        }

        write_field(dst, tag, &value)?;
    }

    Ok(())
}

fn is_valid(src: &[u8]) -> io::Result<bool> {
    let mut buf = src;
    validate(&mut buf)?;
    Ok(true)
}

fn validate(src: &mut &[u8]) -> io::Result<()> {
    while !src.is_empty() {
        split_off_first_chunk::<2>(src).ok_or_else(unexpected_eof)?;
        let ty = src.split_off_first().ok_or_else(unexpected_eof)?;

        match *ty {
            b'A' | b'c' | b'C' => {
                src.split_off_first().ok_or_else(unexpected_eof)?;
            }
            b's' | b'S' => {
                src.split_off(..mem::size_of::<u16>())
                    .ok_or_else(unexpected_eof)?;
            }
            b'i' | b'I' | b'f' => {
                src.split_off(..mem::size_of::<u32>())
                    .ok_or_else(unexpected_eof)?;
            }
            b'Z' => {
                let i = memchr(NUL, src).ok_or_else(unexpected_eof)?;

                // SAFETY: `i < src.len()`.
                let buf = src.split_off(..=i).unwrap();

                // SAFETY: `buf` is nonempty.
                if !buf[..buf.len() - 1]
                    .iter()
                    .all(|b| matches!(b, b' '..=b'~'))
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "invalid string",
                    ));
                }
            }
            b'H' => {
                let i = memchr(NUL, src).ok_or_else(unexpected_eof)?;

                // SAFETY: `i < src.len()`.
                let buf = src.split_off(..=i).unwrap();

                // SAFETY: `buf` is nonempty.
                if !buf.len().is_multiple_of(2)
                    || !buf[..buf.len() - 1]
                        .iter()
                        .all(|b| matches!(b, b'0'..=b'9' | b'A'..=b'F'))
                {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid hex"));
                }
            }
            b'B' => {
                let subtype = src.split_off_first().ok_or_else(unexpected_eof)?;

                let count = split_off_u32_le(src)
                    .ok_or_else(unexpected_eof)
                    .and_then(|n| {
                        usize::try_from(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
                    })?;

                let size = match subtype {
                    b'c' | b'C' => mem::size_of::<u8>(),
                    b's' | b'S' => mem::size_of::<u16>(),
                    b'i' | b'I' | b'f' => mem::size_of::<u32>(),
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
    }

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
    use noodles_sam::alignment::record_buf::{Data as DataBuf, data::field::Value};

    use super::*;

    #[test]
    fn test_write_data() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, data: &DataBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            let d = DataRef::Data(Box::new(data));
            write_data(buf, d)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let data = DataBuf::default();
        t(&mut buf, &data, &[])?;

        let data = [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
            .into_iter()
            .collect();
        t(&mut buf, &data, &[b'N', b'H', b'C', 0x01])?;

        let data = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (Tag::READ_GROUP, Value::from("rg0")),
        ]
        .into_iter()
        .collect();
        t(
            &mut buf,
            &data,
            &[
                b'N', b'H', b'C', 0x01, // NH:C:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_field_encoded_data() -> io::Result<()> {
        let mut dst = Vec::new();

        dst.clear();
        let src = [];
        write_field_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [b'N', b'H', b'C', 0x01]; // NH:C:1
        write_field_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [
            b'N', b'H', b'C', 0x01, // NH:C:1
            b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
        ];
        write_field_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        Ok(())
    }
}
