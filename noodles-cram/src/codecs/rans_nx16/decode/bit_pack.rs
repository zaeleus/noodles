use std::{io, num::NonZero};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub struct Context<'a> {
    pub symbol_count: NonZero<usize>,
    pub mapping_table: &'a [u8],
    pub uncompressed_size: usize,
}

pub fn read_context<'a>(
    src: &mut &'a [u8],
    uncompressed_size: usize,
) -> io::Result<(Context<'a>, usize)> {
    let symbol_count = read_u8(src).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let (mapping_table, rest) = src
        .split_at_checked(symbol_count.get())
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    let len = read_uint7_as(src)?;

    let context = Context {
        symbol_count,
        mapping_table,
        uncompressed_size,
    };

    Ok((context, len))
}

pub fn decode(
    src: &[u8],
    Context {
        symbol_count: n_sym,
        mapping_table: p,
        uncompressed_size: len,
    }: &Context<'_>,
) -> io::Result<Vec<u8>> {
    let mut dst = vec![0; *len];
    let mut j = 0;

    if n_sym.get() <= 1 {
        dst.fill(p[0]);
    } else if n_sym.get() <= 2 {
        let mut v = 0;

        for (i, b) in dst.iter_mut().enumerate() {
            if i % 8 == 0 {
                v = src[j];
                j += 1;
            }

            *b = p[usize::from(v & 0x01)];
            v >>= 1;
        }
    } else if n_sym.get() <= 4 {
        let mut v = 0;

        for (i, b) in dst.iter_mut().enumerate() {
            if i % 4 == 0 {
                v = src[j];
                j += 1;
            }

            *b = p[usize::from(v & 0x03)];
            v >>= 2;
        }
    } else if n_sym.get() <= 16 {
        let mut v = 0;

        for (i, b) in dst.iter_mut().enumerate() {
            if i % 2 == 0 {
                v = src[j];
                j += 1;
            }

            *b = p[usize::from(v & 0x0f)];
            v >>= 4;
        }
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("expected n_sym to be <= 16, got {n_sym}"),
        ));
    }

    Ok(dst)
}
