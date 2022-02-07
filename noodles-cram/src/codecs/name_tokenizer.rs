#![allow(dead_code)]

use std::io::{self, BufRead, Cursor, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use super::rans_nx16::rans_decode_nx16;
use crate::reader::num::read_uint7;

pub fn decode_names<R>(reader: &mut R) -> io::Result<Vec<String>>
where
    R: Read,
{
    let _ulen = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_names = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let use_arith = reader.read_u8()? == 1;

    let mut b = decode_token_byte_streams(reader, use_arith, n_names)?;

    let mut names = vec![String::new(); n_names];
    let mut tokens = vec![vec![None; 128]; n_names];

    for i in 0..n_names {
        decode_single_name(&mut b, &mut names, &mut tokens, i)?;
    }

    Ok(names)
}

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Type {
    Type,
    String,
    Char,
    Digits0,
    DZLen,
    Dup,
    Diff,
    Digits,
    Delta,
    Delta0,
    Match,
    Nop,
    End,
}

impl TryFrom<u8> for Type {
    type Error = String;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n & 0x3f {
            0 => Ok(Self::Type),
            1 => Ok(Self::String),
            2 => Ok(Self::Char),
            3 => Ok(Self::Digits0),
            4 => Ok(Self::DZLen),
            5 => Ok(Self::Dup),
            6 => Ok(Self::Diff),
            7 => Ok(Self::Digits),
            8 => Ok(Self::Delta),
            9 => Ok(Self::Delta0),
            10 => Ok(Self::Match),
            11 => Ok(Self::Nop),
            12 => Ok(Self::End),
            _ => Err(format!("invalid type: expected <= 12, got {}", n)),
        }
    }
}

impl From<Type> for u8 {
    fn from(ty: Type) -> Self {
        match ty {
            Type::Type => 0,
            Type::String => 1,
            Type::Char => 2,
            Type::Digits0 => 3,
            Type::DZLen => 4,
            Type::Dup => 5,
            Type::Diff => 6,
            Type::Digits => 7,
            Type::Delta => 8,
            Type::Delta0 => 9,
            Type::Match => 10,
            Type::Nop => 11,
            Type::End => 12,
        }
    }
}

impl From<Type> for usize {
    fn from(ty: Type) -> Self {
        Self::from(u8::from(ty))
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum Token {
    Char(char),
    String(String),
    Digits(u32),
    PaddedDigits(u32, u8),
    Nop,
}

#[derive(Clone, Debug, Default)]
struct TokenReader {
    type_reader: Cursor<Vec<u8>>,
    string_reader: Cursor<Vec<u8>>,
    char_reader: Cursor<Vec<u8>>,
    digits0_reader: Cursor<Vec<u8>>,
    dz_len_reader: Cursor<Vec<u8>>,
    dup_reader: Cursor<Vec<u8>>,
    diff_reader: Cursor<Vec<u8>>,
    digits_reader: Cursor<Vec<u8>>,
    delta_reader: Cursor<Vec<u8>>,
    delta0_reader: Cursor<Vec<u8>>,
}

impl TokenReader {
    fn get(&self, ty: Type) -> &Cursor<Vec<u8>> {
        match ty {
            Type::Type => &self.type_reader,
            Type::String => &self.string_reader,
            Type::Char => &self.char_reader,
            Type::Digits0 => &self.digits0_reader,
            Type::Dup => &self.dup_reader,
            Type::Diff => &self.diff_reader,
            Type::Digits => &self.digits_reader,
            Type::Delta => &self.delta_reader,
            Type::Delta0 => &self.delta0_reader,
            _ => unimplemented!("unhandled ty: {:?}", ty),
        }
    }

    fn get_mut(&mut self, ty: Type) -> &mut Cursor<Vec<u8>> {
        match ty {
            Type::Type => &mut self.type_reader,
            Type::String => &mut self.string_reader,
            Type::Char => &mut self.char_reader,
            Type::Digits0 => &mut self.digits0_reader,
            Type::Diff => &mut self.diff_reader,
            Type::DZLen => &mut self.dz_len_reader,
            Type::Digits => &mut self.digits_reader,
            Type::Delta => &mut self.delta_reader,
            Type::Delta0 => &mut self.delta0_reader,
            _ => unimplemented!("unhandled ty: {:?}", ty),
        }
    }

    fn set(&mut self, ty: Type, buf: Vec<u8>) {
        match ty {
            Type::Type => *self.type_reader.get_mut() = buf,
            Type::String => *self.string_reader.get_mut() = buf,
            Type::Char => *self.char_reader.get_mut() = buf,
            Type::Digits0 => *self.digits0_reader.get_mut() = buf,
            Type::Dup => *self.dup_reader.get_mut() = buf,
            Type::Diff => *self.diff_reader.get_mut() = buf,
            Type::DZLen => *self.dz_len_reader.get_mut() = buf,
            Type::Digits => *self.digits_reader.get_mut() = buf,
            Type::Delta => *self.delta_reader.get_mut() = buf,
            Type::Delta0 => *self.delta0_reader.get_mut() = buf,
            _ => unimplemented!("unhandled ty: {:?}", ty),
        }
    }

    fn read_type(&mut self) -> io::Result<Type> {
        self.type_reader.read_u8().and_then(|n| {
            Type::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    fn read_distance(&mut self, ty: Type) -> io::Result<usize> {
        assert!(matches!(ty, Type::Dup | Type::Diff));

        self.get_mut(ty).read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    fn read_token(&mut self, prev_token: Option<&Token>) -> io::Result<Option<Token>> {
        match self.read_type()? {
            Type::Char => {
                let c = self.char_reader.read_u8().map(char::from)?;
                Ok(Some(Token::Char(c)))
            }
            Type::String => {
                let mut buf = Vec::new();
                self.string_reader.read_until(0x00, &mut buf)?;
                buf.pop();

                let s = String::from_utf8(buf)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                Ok(Some(Token::String(s)))
            }
            Type::Digits => {
                let d = self.digits_reader.read_u32::<LittleEndian>()?;
                Ok(Some(Token::Digits(d)))
            }
            Type::Digits0 => {
                let d = self.digits0_reader.read_u32::<LittleEndian>()?;
                let l = self.dz_len_reader.read_u8()?;
                Ok(Some(Token::PaddedDigits(d, l)))
            }
            Type::Delta => {
                let delta = self.delta_reader.read_u8().map(u32::from)?;

                match prev_token {
                    Some(Token::Digits(n)) => Ok(Some(Token::Digits(n + delta))),
                    _ => Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invlaid previous token: {:?}", prev_token),
                    )),
                }
            }
            Type::Delta0 => {
                let delta = self.delta0_reader.read_u8().map(u32::from)?;

                match prev_token {
                    Some(Token::PaddedDigits(n, width)) => {
                        Ok(Some(Token::PaddedDigits(n + delta, *width)))
                    }
                    _ => Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invlaid previous token: {:?}", prev_token),
                    )),
                }
            }
            Type::Match => Ok(prev_token.cloned()),
            Type::End => Ok(None),
            _ => Ok(Some(Token::Nop)),
        }
    }
}

fn decode_token_byte_streams<R>(
    reader: &mut R,
    use_arith: bool,
    n_names: usize,
) -> io::Result<Vec<TokenReader>>
where
    R: Read,
{
    let mut b = Vec::new();
    let mut t = -1;

    loop {
        let ttype = match reader.read_u8() {
            Ok(n) => n,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        };

        let tok_new = ttype & 0x80 != 0;
        let tok_dup = ttype & 0x40 != 0;

        let ty =
            Type::try_from(ttype).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if tok_new {
            t += 1;

            b.push(TokenReader::default());

            if ty != Type::Type {
                let mut buf = vec![u8::from(Type::Match); n_names];
                buf[0] = u8::from(ty);
                b[t as usize].set(ty, buf);
            }
        }

        if tok_dup {
            let dup_pos = reader.read_u8().map(usize::from)?;

            let dup_type = reader.read_u8().and_then(|n| {
                Type::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

            let buf = b[dup_pos].get(dup_type).get_ref().clone();
            b[t as usize].set(ty, buf);
        } else {
            let clen = read_uint7(reader).and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

            let mut data = vec![0; clen];
            reader.read_exact(&mut data)?;

            if use_arith {
                todo!("arith_decode");
            } else {
                let mut data_reader = &data[..];
                let buf = rans_decode_nx16(&mut data_reader, 0)?;
                b[t as usize].set(ty, buf);
            }
        }
    }

    Ok(b)
}

fn decode_single_name(
    b: &mut Vec<TokenReader>,
    names: &mut Vec<String>,
    tokens: &mut Vec<Vec<Option<Token>>>,
    n: usize,
) -> io::Result<String> {
    let ty = b[0].read_type()?;
    let dist = b[0].read_distance(ty)?;

    let m = n - dist;

    if ty == Type::Dup {
        names[n] = names[m].clone();
        tokens[n] = tokens[m].clone();
        return Ok(names[n].clone());
    }

    let mut t = 1;

    loop {
        let prev_token = tokens[m][t].as_ref();

        if let Some(token) = b[t].read_token(prev_token)? {
            match &token {
                Token::Char(c) => names[n].push(*c),
                Token::String(s) => names[n].push_str(s),
                Token::Digits(d) => names[n].push_str(&d.to_string()),
                Token::PaddedDigits(d, l) => {
                    let s = format!("{:0width$}", d, width = usize::from(*l));
                    names[n].push_str(&s);
                }
                Token::Nop => {}
            }

            tokens[n][t] = Some(token);
        } else {
            break;
        }

        t += 1;
    }

    Ok(names[n].clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_names() -> io::Result<()> {
        let data = [
            0x58, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x06,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x06, 0x18, 0x00, 0x0c, 0x00, 0x01, 0x00, 0x00, 0x0e, 0x02,
            0x00, 0x22, 0x25, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc,
            0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x01, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02,
            0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x01,
            0x1b, 0x00, 0x04, 0x00, 0x31, 0x37, 0x49, 0x00, 0x01, 0x01, 0x01, 0x01, 0x00, 0x0c,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x5f, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x0a, 0x00, 0x01,
            0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x03, 0x19, 0x00, 0x04, 0x00, 0x22, 0x3d, 0x00, 0x02, 0x01, 0x01,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
            0x01, 0x00, 0x04, 0x15, 0x00, 0x01, 0x05, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00,
            0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x17, 0x00, 0x04, 0x00, 0x02, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00,
            0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x03, 0x01, 0x00, 0xa8,
            0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x1a, 0x00, 0x08, 0x00, 0x7b, 0x7c, 0x00, 0x00, 0x06, 0x01, 0x01, 0x00, 0x7c,
            0x20, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x3a, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x07, 0x00, 0x04, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x07, 0x22, 0x00, 0x0c, 0x00, 0x06, 0x2d, 0x64, 0x65, 0x00, 0xb2, 0xf0, 0x00,
            0x0a, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0xcd, 0x0b, 0x08, 0x00, 0xaf, 0x0e,
            0x08, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x07, 0x00, 0x03, 0x01, 0x00, 0xa8, 0x00, 0x00,
            0x00, 0xa8, 0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x03, 0x1d,
            0x00, 0x08, 0x00, 0x06, 0x21, 0xa3, 0xe3, 0x00, 0x04, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x6e, 0x20, 0x00, 0x00, 0x58, 0x20, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02,
            0x00, 0x04, 0x15, 0x00, 0x02, 0x05, 0x00, 0x02, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x19, 0x00, 0x04,
            0x00, 0x21, 0x3f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x08, 0x02, 0x00, 0x00, 0x0c, 0x02,
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x23, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x17,
            0x00, 0x04, 0x00, 0x09, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00,
            0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x0c,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00,
        ];

        let mut reader = &data[..];
        let actual = decode_names(&mut reader)?;

        let expected = vec![
            String::from("I17_08765:2:123:61541:01763#9"),
            String::from("I17_08765:2:123:1636:08611#9"),
            String::from("I17_08765:2:124:45613:16161#9"),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
