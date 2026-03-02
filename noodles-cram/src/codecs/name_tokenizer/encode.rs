use std::{
    collections::HashMap,
    io::{self, Write},
};

use super::Type;
use crate::io::writer::num::{write_u8, write_u32_le, write_uint7};

const NUL: u8 = 0x00;

/// Compression method for name tokenizer token byte streams.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CompressionMethod {
    RansNx16,
    Aac,
}

pub fn encode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    if let Some(buf) = src.strip_suffix(&[NUL]) {
        src = buf;
    }

    let names: Vec<_> = src.split(|&b| b == NUL).collect();

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

    // Build token writers (one for mode stream, one per token position)
    let mut token_writers = Vec::with_capacity(max_token_count + 1);

    let mut mode_writer = TokenWriter::default();
    for diff in &diffs {
        let token = match diff.mode {
            Mode::Dup(delta) => Token::Dup(delta),
            Mode::Diff(delta) => Token::Diff(delta),
        };
        mode_writer.write_token(&token)?;
    }
    token_writers.push(mode_writer);

    for i in 0..max_token_count {
        let mut tw = TokenWriter::default();
        for diff in &diffs {
            if diff.is_dup() {
                continue;
            }
            if let Some(token) = diff.tokens.get(i) {
                tw.write_token(token)?;
            }
        }
        token_writers.push(tw);
    }

    // Encode all streams with rans, measure total size
    let mut rans_dst = Vec::new();
    write_header(
        &mut rans_dst,
        src.len(),
        names.len(),
        CompressionMethod::RansNx16,
    )?;
    for (pos, tw) in token_writers.iter().enumerate() {
        encode_token_byte_streams(
            &mut rans_dst,
            tw,
            CompressionMethod::RansNx16,
            &token_writers,
            pos,
        )?;
    }

    // Encode all streams with AAC, measure total size
    let mut aac_dst = Vec::new();
    write_header(&mut aac_dst, src.len(), names.len(), CompressionMethod::Aac)?;
    for (pos, tw) in token_writers.iter().enumerate() {
        encode_token_byte_streams(
            &mut aac_dst,
            tw,
            CompressionMethod::Aac,
            &token_writers,
            pos,
        )?;
    }

    // Pick whichever is smaller
    if aac_dst.len() < rans_dst.len() {
        Ok(aac_dst)
    } else {
        Ok(rans_dst)
    }
}

fn write_header<W>(
    writer: &mut W,
    src_len: usize,
    names_count: usize,
    method: CompressionMethod,
) -> io::Result<()>
where
    W: Write,
{
    let ulen =
        u32::try_from(src_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, ulen)?;

    let n_names =
        u32::try_from(names_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n_names)?;

    let use_arith: u8 = match method {
        CompressionMethod::RansNx16 => 0,
        CompressionMethod::Aac => 1,
    };
    write_u8(writer, use_arith)?;

    Ok(())
}

#[derive(Debug)]
enum Token {
    String(Vec<u8>),
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
    raw_tokens: Vec<Vec<u8>>,
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

fn tokenize(b: &[u8]) -> impl Iterator<Item = &[u8]> {
    use std::iter;

    let mut start = 0;
    let mut end = 0;

    iter::from_fn(move || {
        while end < b.len() && b[end].is_ascii_alphanumeric() {
            end += 1;
        }

        if start != end {
            let beg = start;
            start = end;
            return Some(&b[beg..end]);
        }

        while end < b.len() && !b[end].is_ascii_alphanumeric() {
            end += 1;
        }

        if start == end {
            None
        } else {
            let beg = start;
            start = end;
            Some(&b[beg..end])
        }
    })
}

fn build_first_diff(name: &[u8]) -> Diff {
    let mut diff = Diff::new(Mode::Diff(0));

    let raw_tokens = tokenize(name);

    for raw_token in raw_tokens {
        let token = if let Some(n) = parse_digits0(raw_token) {
            Token::PaddedDigits(n, raw_token.len())
        } else if let Some(n) = parse_digits(raw_token) {
            Token::Digits(n)
        } else if raw_token.len() == 1 {
            Token::Char(raw_token[0])
        } else {
            Token::String(raw_token.into())
        };

        diff.raw_tokens.push(raw_token.into());
        diff.tokens.push(token);
    }

    diff.tokens.push(Token::End);

    diff
}

fn build_diff(
    diffs: &[Diff],
    names_indices: &HashMap<&[u8], usize>,
    i: usize,
    name: &[u8],
) -> Diff {
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
                Some(Token::Char(raw_token[0]))
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

fn parse_digits0(s: &[u8]) -> Option<u32> {
    if s.starts_with(b"0") {
        parse_u32(s).ok()
    } else {
        None
    }
}

fn parse_digits(s: &[u8]) -> Option<u32> {
    parse_u32(s).ok()
}

fn parse_delta(prev_token: &Token, s: &[u8]) -> Option<(u32, u8)> {
    if let Token::Digits(n) | Token::Delta(n, _) = prev_token {
        let m = parse_u32(s).ok()?;

        if m >= *n {
            let delta = m - n;

            if delta <= u32::from(u8::MAX) {
                return Some((m, delta as u8));
            }
        }
    }

    None
}

fn parse_delta0(prev_s: &[u8], prev_token: &Token, s: &[u8]) -> Option<(u32, u8)> {
    if let Token::PaddedDigits(n, _) | Token::Delta0(n, _) = prev_token
        && s.len() == prev_s.len()
    {
        let m = parse_u32(s).ok()?;

        if m >= *n {
            let delta = m - n;

            if delta <= u32::from(u8::MAX) {
                return Some((m, delta as u8));
            }
        }
    }

    None
}

fn parse_u32(src: &[u8]) -> io::Result<u32> {
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
                self.string_writer.write_all(s)?;
                self.string_writer.write_all(&[NUL])?;
            }
            Token::Char(b) => write_u8(&mut self.char_writer, *b)?,
            Token::PaddedDigits(n, width) => {
                // d
                write_u32_le(&mut self.digits0_writer, *n)?;

                let l = u8::try_from(*width)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u8(&mut self.dz_len_writer, l)?;
            }
            // ...
            Token::Dup(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u32_le(&mut self.dup_writer, n)?;
            }
            Token::Diff(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u32_le(&mut self.diff_writer, n)?;
            }
            Token::Digits(n) => {
                write_u32_le(&mut self.digits_writer, *n)?;
            }
            Token::Delta(_, delta) => {
                write_u8(&mut self.delta_writer, *delta)?;
            }
            Token::Delta0(_, delta) => {
                write_u8(&mut self.delta0_writer, *delta)?;
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

    write_u8(writer, n)
}

/// Returns all (type, data) pairs for a given TokenWriter in canonical order.
fn token_writer_streams(tw: &TokenWriter) -> [(Type, &[u8]); 10] {
    [
        (Type::Type, tw.type_writer.as_slice()),
        (Type::String, tw.string_writer.as_slice()),
        (Type::Char, tw.char_writer.as_slice()),
        (Type::Digits0, tw.digits0_writer.as_slice()),
        (Type::DZLen, tw.dz_len_writer.as_slice()),
        (Type::Dup, tw.dup_writer.as_slice()),
        (Type::Diff, tw.diff_writer.as_slice()),
        (Type::Digits, tw.digits_writer.as_slice()),
        (Type::Delta, tw.delta_writer.as_slice()),
        (Type::Delta0, tw.delta0_writer.as_slice()),
    ]
}

fn encode_token_byte_streams<W>(
    writer: &mut W,
    token_writer: &TokenWriter,
    method: CompressionMethod,
    all_token_writers: &[TokenWriter],
    current_pos: usize,
) -> io::Result<()>
where
    W: Write,
{
    for (ty, buf) in token_writer_streams(token_writer) {
        if buf.is_empty() {
            continue;
        }

        // Check for tok_dup: see if an earlier position has an identical byte stream
        if let Some((dup_pos, dup_type)) =
            find_duplicate_stream(all_token_writers, current_pos, buf)
        {
            // Emit type byte with 0x40 (tok_dup) flag set
            let type_byte = if ty == Type::Type {
                0x80 | 0x40
            } else {
                u8::from(ty) | 0x40
            };
            write_u8(writer, type_byte)?;
            let dup_pos_u8 = u8::try_from(dup_pos)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            write_u8(writer, dup_pos_u8)?;
            write_u8(writer, u8::from(dup_type))?;
        } else {
            encode_token_byte_stream(writer, ty, buf, method)?;
        }
    }

    Ok(())
}

/// Search all earlier token positions for a byte stream identical to `buf`.
fn find_duplicate_stream(
    all_token_writers: &[TokenWriter],
    current_pos: usize,
    buf: &[u8],
) -> Option<(usize, Type)> {
    for (pos, tw) in all_token_writers[..current_pos].iter().enumerate() {
        for (ty, stream) in token_writer_streams(tw) {
            if !stream.is_empty() && stream == buf {
                return Some((pos, ty));
            }
        }
    }
    None
}

fn encode_token_byte_stream<W>(
    writer: &mut W,
    ty: Type,
    buf: &[u8],
    method: CompressionMethod,
) -> io::Result<()>
where
    W: Write,
{
    use crate::codecs::{aac, rans_nx16};

    if ty == Type::Type {
        write_u8(writer, 0x80)?;
    } else {
        write_type(writer, ty)?;
    }

    let cdata = match method {
        CompressionMethod::RansNx16 => rans_nx16::encode(rans_nx16::Flags::empty(), buf)?,
        CompressionMethod::Aac => aac::encode(aac::Flags::empty(), buf)?,
    };

    let clen =
        u32::try_from(cdata.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(writer, clen)?;

    writer.write_all(&cdata)?;

    Ok(())
}
