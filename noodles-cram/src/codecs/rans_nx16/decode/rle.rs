use std::{borrow::Cow, io};

use super::{order_0, split_off};
use crate::io::reader::num::{read_u8, read_uint7, read_uint7_as};

pub(super) struct Context<'a> {
    src: Cow<'a, [u8]>,
    len: usize,
}

pub(super) fn read_context<'a>(
    src: &mut &'a [u8],
    state_count: usize,
    uncompressed_size: usize,
) -> io::Result<(Context<'a>, usize)> {
    let (context_size, is_compressed) = read_header(src)?;

    let len = read_uint7_as(src)?;
    let context_src = read_src(src, state_count, context_size, is_compressed)?;

    let ctx = Context {
        src: context_src,
        len: uncompressed_size,
    };

    Ok((ctx, len))
}

fn read_header(src: &mut &[u8]) -> io::Result<(usize, bool)> {
    let n = read_uint7(src)?;

    let context_size =
        usize::try_from(n >> 1).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let is_compressed = (n & 0x01) == 0;

    Ok((context_size, is_compressed))
}

fn read_src<'a>(
    src: &mut &'a [u8],
    state_count: usize,
    len: usize,
    is_compressed: bool,
) -> io::Result<Cow<'a, [u8]>> {
    if is_compressed {
        let compressed_size = read_uint7_as(src)?;
        let mut buf = split_off(src, compressed_size)?;

        let mut dst = vec![0; len];
        order_0::decode(&mut buf, &mut dst, state_count)?;
        Ok(Cow::from(dst))
    } else {
        split_off(src, len).map(Cow::from)
    }
}

pub fn decode(mut src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    let mut context_src = &ctx.src[..];

    let l = read_rle_table(&mut context_src)?;

    let mut dst = vec![0; ctx.len];
    let mut j = 0;

    while j < dst.len() {
        let sym = read_u8(&mut src)?;

        if l[usize::from(sym)] {
            let run = read_uint7_as(&mut context_src)?;

            for k in 0..=run {
                dst[j + k] = sym;
            }

            j += run + 1;
        } else {
            dst[j] = sym;
            j += 1;
        }
    }

    Ok(dst)
}

fn read_rle_table(src: &mut &[u8]) -> io::Result<[bool; 256]> {
    let mut rle_table = [false; 256];

    let symbol_count = {
        let n = read_u8(src).map(usize::from)?;
        if n == 0 { rle_table.len() } else { n }
    };

    for _ in 0..symbol_count {
        let sym = read_u8(src)?;
        rle_table[usize::from(sym)] = true;
    }

    Ok(rle_table)
}
