mod header;

use std::{
    array,
    io::{self, Write},
    num::NonZero,
};

use bstr::ByteSlice;

use self::header::{CompressionMethod, read_header};
use super::{TYPE_COUNT, Type};
use crate::{
    codecs::{aac, rans_nx16},
    io::reader::num::{read_u8, read_u32_le, read_uint7_as},
};

// § 5 "Name tokenisation codec" (2025-04-07): "A maximum of 128 tokens are permitted within any
// single read name."
const MAX_TOKEN_COUNT: usize = 128;

const NUL: u8 = 0x00;

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let (uncompressed_size, name_count, compression_method) = read_header(&mut src)?;

    let bufs = read_buffers(&mut src, name_count, compression_method)?;
    let mut readers = build_readers(&bufs);

    let mut tokens = Vec::with_capacity(name_count.get());

    for i in 0..name_count.get() {
        let ts = read_tokens(&mut readers, &tokens, i)?;
        tokens.push(ts);
    }

    let mut dst = Vec::with_capacity(uncompressed_size);
    write_names(&mut dst, &tokens)?;

    Ok(dst)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Token<'a> {
    Char(u8),
    String(&'a [u8]),
    Digits(u32),
    PaddedDigits(u32, u8),
    Nop,
}

#[derive(Clone, Debug, Default)]
struct TokenReader<'a>([&'a [u8]; TYPE_COUNT]);

impl<'a> TokenReader<'a> {
    fn read_type(&mut self) -> io::Result<Type> {
        read_u8(&mut self.0[usize::from(Type::Type)]).and_then(decode_type)
    }

    fn read_distance(&mut self, ty: Type) -> io::Result<usize> {
        let reader = match ty {
            Type::Dup => &mut self.0[usize::from(Type::Dup)],
            Type::Diff => &mut self.0[usize::from(Type::Diff)],
            _ => return Err(invalid_type(ty)),
        };

        read_u32_le(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    fn read_token(&mut self, prev_token: Option<Token<'a>>) -> io::Result<Option<Token<'a>>> {
        match self.read_type()? {
            Type::Char => read_char(&mut self.0[usize::from(Type::Char)]).map(Some),
            Type::String => read_string(&mut self.0[usize::from(Type::String)]).map(Some),
            Type::Digits => read_digits(&mut self.0[usize::from(Type::Digits)]).map(Some),
            Type::Digits0 => {
                let n = read_u32_le(&mut self.0[usize::from(Type::Digits0)])?;
                let width = read_u8(&mut self.0[usize::from(Type::DZLen)])?;
                Ok(Some(Token::PaddedDigits(n, width)))
            }
            Type::Delta => read_delta(&mut self.0[usize::from(Type::Delta)], prev_token).map(Some),
            Type::Delta0 => {
                read_delta0(&mut self.0[usize::from(Type::Delta0)], prev_token).map(Some)
            }
            Type::Match => prev_token
                .filter(|&token| token != Token::Nop)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid match: invalid previous token: {prev_token:?}"),
                    )
                })
                .map(Some),
            Type::Nop => Ok(Some(Token::Nop)),
            Type::End => Ok(None),
            ty => Err(invalid_type(ty)),
        }
    }
}

fn read_string<'a>(src: &mut &'a [u8]) -> io::Result<Token<'a>> {
    let i = src
        .as_bstr()
        .find_byte(NUL)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let (buf, rest) = src.split_at(i);

    *src = &rest[1..];

    Ok(Token::String(buf))
}

fn read_char<'a>(src: &mut &'a [u8]) -> io::Result<Token<'a>> {
    read_u8(src).map(Token::Char)
}

fn read_digits<'a>(src: &mut &'a [u8]) -> io::Result<Token<'a>> {
    read_u32_le(src).map(Token::Digits)
}

fn read_delta<'a>(src: &mut &'a [u8], prev_token: Option<Token<'a>>) -> io::Result<Token<'a>> {
    let delta = read_u8(src).map(u32::from)?;

    if let Some(Token::Digits(n)) = prev_token {
        n.checked_add(delta).map(Token::Digits).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid delta: attempt to add with overflow",
            )
        })
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid previous token: {prev_token:?}"),
        ))
    }
}

fn read_delta0<'a>(src: &mut &'a [u8], prev_token: Option<Token<'a>>) -> io::Result<Token<'a>> {
    let delta = read_u8(src).map(u32::from)?;

    if let Some(Token::PaddedDigits(n, width)) = prev_token {
        n.checked_add(delta)
            .map(|m| Token::PaddedDigits(m, width))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid delta0: attempt to add with overflow",
                )
            })
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid previous token: {prev_token:?}"),
        ))
    }
}

#[derive(Default)]
enum Buffer {
    #[default]
    Empty,
    Data(Vec<u8>),
    Duplicate(usize, Type),
}

fn read_buffers(
    src: &mut &[u8],
    name_count: NonZero<usize>,
    compression_method: CompressionMethod,
) -> io::Result<Vec<[Buffer; TYPE_COUNT]>> {
    let mut bufs = Vec::new();

    while !src.is_empty() {
        let n = read_u8(src)?;

        let is_new = n & 0x80 != 0;
        let is_dup = n & 0x40 != 0;
        let ty = decode_type(n & 0x3f)?;

        if is_new {
            if bufs.len() >= MAX_TOKEN_COUNT {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid buffer count",
                ));
            }

            let mut buf = array::from_fn(|_| Buffer::default());

            if ty != Type::Type {
                let mut b = vec![u8::from(Type::Match); name_count.get()];
                b[0] = u8::from(ty);
                buf[0] = Buffer::Data(b);
            }

            bufs.push(buf);
        }

        let buf = bufs
            .last()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty buffers"))?;

        if !matches!(buf[usize::from(ty)], Buffer::Empty) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate buffer",
            ));
        }

        let ty_buf = if is_dup {
            let j = read_u8(src).map(usize::from)?;
            let dup_type = read_u8(src).and_then(decode_type)?;

            let buf = bufs.get(j).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid duplicate buffer position",
                )
            })?;

            match buf[usize::from(dup_type)] {
                Buffer::Empty => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid duplicate buffer: empty source",
                    ));
                }
                Buffer::Data(_) => Buffer::Duplicate(j, dup_type),
                Buffer::Duplicate(k, orig_type) => Buffer::Duplicate(k, orig_type),
            }
        } else {
            let compressed_size = read_uint7_as(src)?;
            let buf = split_off(src, compressed_size)?;

            let uncompressed_data = match compression_method {
                CompressionMethod::RansNx16 => rans_nx16::decode(buf, 0)?,
                CompressionMethod::AdaptiveArithmeticCoder => aac::decode(buf, 0)?,
            };

            Buffer::Data(uncompressed_data)
        };

        let buf = bufs
            .last_mut()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty buffers"))?;

        buf[usize::from(ty)] = ty_buf;
    }

    Ok(bufs)
}

fn decode_type(n: u8) -> io::Result<Type> {
    match n {
        0 => Ok(Type::Type),
        1 => Ok(Type::String),
        2 => Ok(Type::Char),
        3 => Ok(Type::Digits0),
        4 => Ok(Type::DZLen),
        5 => Ok(Type::Dup),
        6 => Ok(Type::Diff),
        7 => Ok(Type::Digits),
        8 => Ok(Type::Delta),
        9 => Ok(Type::Delta0),
        10 => Ok(Type::Match),
        11 => Ok(Type::Nop),
        12 => Ok(Type::End),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid type: expected <= 12, got {n}"),
        )),
    }
}

fn build_readers(bufs: &[[Buffer; TYPE_COUNT]]) -> Vec<TokenReader<'_>> {
    bufs.iter()
        .map(|bs| {
            bs.each_ref().map(|buf| match buf {
                Buffer::Empty => &[][..],
                Buffer::Data(b) => b.as_slice(),
                Buffer::Duplicate(j, ty) => match &bufs[*j][usize::from(*ty)] {
                    Buffer::Data(b) => b.as_slice(),
                    _ => unreachable!(),
                },
            })
        })
        .map(TokenReader)
        .collect()
}

fn read_tokens<'a>(
    readers: &mut [TokenReader<'a>],
    tokens: &[Vec<Token<'a>>],
    i: usize,
) -> io::Result<Vec<Token<'a>>> {
    let (meta_reader, token_readers) = readers
        .split_first_mut()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty readers"))?;

    let ty = meta_reader.read_type()?;

    let dist = meta_reader.read_distance(ty)?;

    let j = i.checked_sub(dist).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid distance: attempt to subtract with overflow",
        )
    })?;

    if ty == Type::Dup {
        tokens
            .get(j)
            .cloned()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid duplicate distance"))
    } else {
        let prev_tokens = tokens.get(j);
        let mut ts = Vec::with_capacity(token_readers.len());

        for (t, reader) in token_readers.iter_mut().enumerate() {
            let prev_token = prev_tokens.and_then(|ts| ts.get(t).copied());

            if let Some(token) = reader.read_token(prev_token)? {
                ts.push(token);
            } else {
                return Ok(ts);
            }
        }

        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing end token",
        ))
    }
}

fn write_names(dst: &mut Vec<u8>, tokens: &[Vec<Token>]) -> io::Result<()> {
    for ts in tokens {
        write_name(dst, ts)?;
        dst.push(NUL);
    }

    Ok(())
}

fn write_name(dst: &mut Vec<u8>, tokens: &[Token]) -> io::Result<()> {
    for token in tokens {
        match token {
            Token::Char(c) => dst.push(*c),
            Token::String(s) => dst.extend(*s),
            Token::Digits(d) => write!(dst, "{d}")?,
            Token::PaddedDigits(d, l) => write!(dst, "{:0width$}", d, width = usize::from(*l))?,
            Token::Nop => {}
        }
    }

    Ok(())
}

fn split_off<'a>(src: &mut &'a [u8], len: usize) -> io::Result<&'a [u8]> {
    src.split_off(..len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn invalid_type(ty: Type) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("invalid type: {ty:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        let src = [
            0x59, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x06, 0x00, 0x06, 0x08, 0xa0, 0x0c, 0x02, 0x00, 0x01, 0x02, 0x10, 0x01, 0x81, 0x0d,
            0x20, 0x0b, 0x49, 0x31, 0x37, 0x5f, 0x30, 0x38, 0x37, 0x36, 0x35, 0x3a, 0x00, 0x87,
            0x06, 0x20, 0x04, 0x02, 0x00, 0x00, 0x00, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05,
            0x20, 0x03, 0x07, 0x0a, 0x08, 0x07, 0x06, 0x20, 0x04, 0x7b, 0x00, 0x00, 0x00, 0x08,
            0x03, 0x20, 0x01, 0x01, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x07, 0x00, 0x07, 0x0e, 0x20, 0x0c, 0x65, 0xf0, 0x00, 0x00, 0x64, 0x06, 0x00, 0x00,
            0x2d, 0xb2, 0x00, 0x00, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x03, 0x00, 0x03, 0x0e, 0x20, 0x0c, 0xe3, 0x06, 0x00, 0x00, 0xa3, 0x21, 0x00, 0x00,
            0x21, 0x3f, 0x00, 0x00, 0x04, 0x05, 0x20, 0x03, 0x05, 0x05, 0x05, 0x82, 0x03, 0x20,
            0x01, 0x23, 0x87, 0x06, 0x20, 0x04, 0x09, 0x00, 0x00, 0x00, 0x80, 0x05, 0xa0, 0x03,
            0x01, 0x0c, 0x00,
        ];

        let actual = decode(&src)?;

        let expected = b"\
I17_08765:2:123:61541:01763#9\0\
I17_08765:2:123:1636:08611#9\0\
I17_08765:2:124:45613:16161#9\0\
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_with_arithmetic_coder() -> io::Result<()> {
        let src = [
            0x59, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x01, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x06, 0x00, 0x06, 0x08, 0xe0, 0x0c, 0x02, 0x00, 0x01, 0x02, 0x10, 0x01, 0x81, 0x0d,
            0x20, 0x0b, 0x49, 0x31, 0x37, 0x5f, 0x30, 0x38, 0x37, 0x36, 0x35, 0x3a, 0x00, 0x87,
            0x06, 0x20, 0x04, 0x02, 0x00, 0x00, 0x00, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05,
            0x20, 0x03, 0x07, 0x0a, 0x08, 0x07, 0x06, 0x20, 0x04, 0x7b, 0x00, 0x00, 0x00, 0x08,
            0x03, 0x20, 0x01, 0x01, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x07, 0x00, 0x07, 0x0e, 0x20, 0x0c, 0x65, 0xf0, 0x00, 0x00, 0x64, 0x06, 0x00, 0x00,
            0x2d, 0xb2, 0x00, 0x00, 0x82, 0x03, 0x20, 0x01, 0x3a, 0x80, 0x05, 0xa0, 0x03, 0x01,
            0x03, 0x00, 0x03, 0x0e, 0x20, 0x0c, 0xe3, 0x06, 0x00, 0x00, 0xa3, 0x21, 0x00, 0x00,
            0x21, 0x3f, 0x00, 0x00, 0x04, 0x05, 0x20, 0x03, 0x05, 0x05, 0x05, 0x82, 0x03, 0x20,
            0x01, 0x23, 0x87, 0x06, 0x20, 0x04, 0x09, 0x00, 0x00, 0x00, 0x80, 0x05, 0xa0, 0x03,
            0x01, 0x0c, 0x00,
        ];

        let actual = decode(&src)?;

        let expected = b"\
I17_08765:2:123:61541:01763#9\0\
I17_08765:2:123:1636:08611#9\0\
I17_08765:2:124:45613:16161#9\0\
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_with_dzlen() -> io::Result<()> {
        let src = [
            0x0a, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x80, 0x04, 0x20, 0x02, 0x06,
            0x06, 0x06, 0x07, 0xa0, 0x08, 0x02, 0x00, 0x01, 0x01, 0x10, 0x80, 0x04, 0x20, 0x02,
            0x03, 0x03, 0x03, 0x05, 0xa0, 0x08, 0x01, 0x00, 0x00, 0x04, 0x04, 0x20, 0x02, 0x01,
            0x02, 0x82, 0x03, 0x20, 0x01, 0x3a, 0xc0, 0x01, 0x00, 0x43, 0x01, 0x03, 0x44, 0x01,
            0x04, 0x80, 0x04, 0x20, 0x02, 0x0c, 0x0c,
        ];

        let actual = decode(&src)?;

        let expected = b"0:0\x0000:00\x00";
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_type() -> io::Result<()> {
        assert_eq!(decode_type(0)?, Type::Type);
        assert_eq!(decode_type(1)?, Type::String);
        assert_eq!(decode_type(2)?, Type::Char);
        assert_eq!(decode_type(3)?, Type::Digits0);
        assert_eq!(decode_type(4)?, Type::DZLen);
        assert_eq!(decode_type(5)?, Type::Dup);
        assert_eq!(decode_type(6)?, Type::Diff);
        assert_eq!(decode_type(7)?, Type::Digits);
        assert_eq!(decode_type(8)?, Type::Delta);
        assert_eq!(decode_type(9)?, Type::Delta0);
        assert_eq!(decode_type(10)?, Type::Match);
        assert_eq!(decode_type(11)?, Type::Nop);
        assert_eq!(decode_type(12)?, Type::End);

        assert!(matches!(
            decode_type(13),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
