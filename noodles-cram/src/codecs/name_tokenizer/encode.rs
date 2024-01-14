use std::{
    collections::HashMap,
    io::{self, Write},
    str,
};

use byteorder::{LittleEndian, WriteBytesExt};

use super::Type;
use crate::io::writer::num::write_uint7;

const NUL: u8 = 0x00;

pub fn encode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    if let Some(buf) = src.strip_suffix(&[NUL]) {
        src = buf;
    }

    let names: Vec<_> = src
        .split(|&b| b == NUL)
        .map(str::from_utf8)
        .collect::<Result<_, _>>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_header(&mut dst, src.len(), names.len())?;

    let mut names_indices = HashMap::new();
    let mut diffs = Vec::with_capacity(names.len());
    let mut max_token_count = 0;

    if let Some(name) = names.first() {
        let diff = build_first_diff(name);
        max_token_count = max_token_count.max(diff.tokens.len());
        diffs.push(diff);
    }

    for (i, name) in names.iter().enumerate().skip(1) {
        let diff = build_diff(&diffs, &names_indices, i, name);
        names_indices.entry(name).or_insert(i);
        max_token_count = max_token_count.max(diff.tokens.len());
        diffs.push(diff);
    }

    let mut token_writer = TokenWriter::default();

    for diff in &diffs {
        let token = match diff.mode {
            Mode::Dup(delta) => Token::Dup(delta),
            Mode::Diff(delta) => Token::Diff(delta),
        };

        token_writer.write_token(&token)?;
    }

    encode_token_byte_streams(&mut dst, &token_writer)?;

    for i in 0..max_token_count {
        let mut token_writer = TokenWriter::default();

        for diff in &diffs {
            if diff.is_dup() {
                continue;
            }

            if let Some(token) = diff.tokens.get(i) {
                token_writer.write_token(token)?;
            }
        }

        encode_token_byte_streams(&mut dst, &token_writer)?;
    }

    Ok(dst)
}

fn write_header<W>(writer: &mut W, src_len: usize, names_count: usize) -> io::Result<()>
where
    W: Write,
{
    let ulen =
        u32::try_from(src_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(ulen)?;

    let n_names =
        u32::try_from(names_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_names)?;

    // TODO: use_arith
    writer.write_u8(0)?;

    Ok(())
}

#[derive(Debug)]
enum Token {
    String(String),
    Char(u8),
    PaddedDigits(u32, usize),
    // ...
    Dup(usize),
    Diff(usize),
    Digits(u32),
    Delta(u32, u8),
    Delta0(u32, u8),
    Match,
    // ...
    End,
}

impl Token {
    fn ty(&self) -> Type {
        match self {
            Self::String(_) => Type::String,
            Self::Char(_) => Type::Char,
            Self::PaddedDigits(..) => Type::Digits0,
            // ...
            Self::Dup(_) => Type::Dup,
            Self::Diff(_) => Type::Diff,
            Self::Digits(_) => Type::Digits,
            Self::Delta(..) => Type::Delta,
            Self::Delta0(..) => Type::Delta0,
            Self::Match => Type::Match,
            // ...
            Self::End => Type::End,
        }
    }
}

enum Mode {
    Diff(usize),
    Dup(usize),
}

struct Diff {
    mode: Mode,
    raw_tokens: Vec<String>,
    tokens: Vec<Token>,
}

impl Diff {
    fn new(mode: Mode) -> Self {
        Self {
            mode,
            raw_tokens: Vec::new(),
            tokens: Vec::new(),
        }
    }

    fn delta(&self) -> usize {
        match self.mode {
            Mode::Diff(delta) => delta,
            Mode::Dup(delta) => delta,
        }
    }

    fn is_dup(&self) -> bool {
        matches!(self.mode, Mode::Dup(_))
    }
}

fn tokenize(s: &str) -> impl Iterator<Item = &str> {
    use std::iter;

    let b = s.as_bytes();
    let mut start = 0;
    let mut end = 0;

    iter::from_fn(move || {
        while end < s.len() && b[end].is_ascii_alphanumeric() {
            end += 1;
        }

        if start != end {
            let beg = start;
            start = end;
            return Some(&s[beg..end]);
        }

        while end < s.len() && !b[end].is_ascii_alphanumeric() {
            end += 1;
        }

        if start == end {
            None
        } else {
            let beg = start;
            start = end;
            Some(&s[beg..end])
        }
    })
}

fn build_first_diff(name: &str) -> Diff {
    let mut diff = Diff::new(Mode::Diff(0));

    let raw_tokens = tokenize(name);

    for raw_token in raw_tokens {
        let token = if let Some(n) = parse_digits0(raw_token) {
            Token::PaddedDigits(n, raw_token.len())
        } else if let Some(n) = parse_digits(raw_token) {
            Token::Digits(n)
        } else if raw_token.len() == 1 {
            let b = raw_token.as_bytes()[0];
            Token::Char(b)
        } else {
            Token::String(raw_token.into())
        };

        diff.raw_tokens.push(raw_token.into());
        diff.tokens.push(token);
    }

    diff.tokens.push(Token::End);

    diff
}

fn build_diff(diffs: &[Diff], names_indices: &HashMap<&str, usize>, i: usize, name: &str) -> Diff {
    let mut diff = if let Some(j) = names_indices.get(name) {
        let delta = i - j;
        Diff::new(Mode::Dup(delta))
    } else {
        Diff::new(Mode::Diff(1))
    };

    let prev_diff = &diffs[i - diff.delta()];
    let raw_tokens = tokenize(name);

    for (j, raw_token) in raw_tokens.enumerate() {
        let mut token = None;

        if let (Some(prev_raw_token), Some(prev_token)) =
            (prev_diff.raw_tokens.get(j), prev_diff.tokens.get(j))
        {
            if raw_token == prev_raw_token {
                token = Some(Token::Match);
            } else if let Some((n, delta)) = parse_delta(prev_token, raw_token) {
                token = Some(Token::Delta(n, delta));
            } else if let Some((n, delta)) = parse_delta0(prev_raw_token, prev_token, raw_token) {
                token = Some(Token::Delta0(n, delta));
            }
        }

        if token.is_none() {
            token = if let Some(n) = parse_digits0(raw_token) {
                Some(Token::PaddedDigits(n, raw_token.len()))
            } else if let Some(n) = parse_digits(raw_token) {
                Some(Token::Digits(n))
            } else if raw_token.len() == 1 {
                let b = raw_token.as_bytes()[0];
                Some(Token::Char(b))
            } else {
                Some(Token::String(raw_token.into()))
            }
        }

        let token = token.unwrap();

        diff.raw_tokens.push(raw_token.into());
        diff.tokens.push(token);
    }

    diff.tokens.push(Token::End);

    diff
}

fn parse_digits0(s: &str) -> Option<u32> {
    if s.starts_with('0') {
        s.parse().ok()
    } else {
        None
    }
}

fn parse_digits(s: &str) -> Option<u32> {
    s.parse().ok()
}

fn parse_delta(prev_token: &Token, s: &str) -> Option<(u32, u8)> {
    if let Token::Digits(n) | Token::Delta(n, _) = prev_token {
        let m = s.parse().ok()?;

        if m >= *n {
            let delta = m - n;

            if delta <= u32::from(u8::MAX) {
                return Some((m, delta as u8));
            }
        }
    }

    None
}

fn parse_delta0(prev_s: &str, prev_token: &Token, s: &str) -> Option<(u32, u8)> {
    if let Token::PaddedDigits(n, _) | Token::Delta0(n, _) = prev_token {
        if s.len() == prev_s.len() {
            let m = s.parse().ok()?;

            if m >= *n {
                let delta = m - n;

                if delta <= u32::from(u8::MAX) {
                    return Some((m, delta as u8));
                }
            }
        }
    }

    None
}

#[derive(Default)]
struct TokenWriter {
    type_writer: Vec<u8>,
    string_writer: Vec<u8>,
    char_writer: Vec<u8>,
    digits0_writer: Vec<u8>,
    dz_len_writer: Vec<u8>,
    dup_writer: Vec<u8>,
    diff_writer: Vec<u8>,
    digits_writer: Vec<u8>,
    delta_writer: Vec<u8>,
    delta0_writer: Vec<u8>,
}

impl TokenWriter {
    fn write_token(&mut self, token: &Token) -> io::Result<()> {
        write_type(&mut self.type_writer, token.ty())?;

        match token {
            Token::String(s) => {
                self.string_writer.write_all(s.as_bytes())?;
                self.string_writer.write_all(&[NUL])?;
            }
            Token::Char(b) => self.char_writer.write_u8(*b)?,
            Token::PaddedDigits(n, width) => {
                // d
                self.digits0_writer.write_u32::<LittleEndian>(*n)?;

                let l = u8::try_from(*width)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                self.dz_len_writer.write_u8(l)?;
            }
            // ...
            Token::Dup(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                self.dup_writer.write_u32::<LittleEndian>(n)?;
            }
            Token::Diff(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                self.diff_writer.write_u32::<LittleEndian>(n)?;
            }
            Token::Digits(n) => {
                self.digits_writer.write_u32::<LittleEndian>(*n)?;
            }
            Token::Delta(_, delta) => {
                self.delta_writer.write_u8(*delta)?;
            }
            Token::Delta0(_, delta) => {
                self.delta0_writer.write_u8(*delta)?;
            }
            Token::Match => {}
            // ...
            Token::End => {}
        }

        Ok(())
    }
}

fn write_type<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    let n = match ty {
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
    };

    writer.write_u8(n)
}

fn encode_token_byte_streams<W>(writer: &mut W, token_writer: &TokenWriter) -> io::Result<()>
where
    W: Write,
{
    encode_token_byte_stream(writer, Type::Type, &token_writer.type_writer)?;
    encode_token_byte_stream(writer, Type::String, &token_writer.string_writer)?;
    encode_token_byte_stream(writer, Type::Char, &token_writer.char_writer)?;
    encode_token_byte_stream(writer, Type::Digits0, &token_writer.digits0_writer)?;
    encode_token_byte_stream(writer, Type::DZLen, &token_writer.dz_len_writer)?;
    encode_token_byte_stream(writer, Type::Dup, &token_writer.dup_writer)?;
    encode_token_byte_stream(writer, Type::Diff, &token_writer.diff_writer)?;
    encode_token_byte_stream(writer, Type::Digits, &token_writer.digits_writer)?;
    encode_token_byte_stream(writer, Type::Delta, &token_writer.delta_writer)?;
    encode_token_byte_stream(writer, Type::Delta0, &token_writer.delta0_writer)?;

    Ok(())
}

fn encode_token_byte_stream<W>(writer: &mut W, ty: Type, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    use crate::codecs::rans_nx16::{self, Flags};

    if buf.is_empty() {
        return Ok(());
    }

    if ty == Type::Type {
        writer.write_u8(0x80)?;
    } else {
        write_type(writer, ty)?;
    }

    let cdata = rans_nx16::encode(Flags::empty(), buf)?;

    let clen =
        u32::try_from(cdata.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(writer, clen)?;

    writer.write_all(&cdata)?;

    Ok(())
}
