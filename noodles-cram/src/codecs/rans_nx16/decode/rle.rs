use std::io::{self, Read};

use super::order_0;
use crate::io::reader::num::{read_u8, read_uint7, read_uint7_as};

pub(super) fn decode_rle_meta(src: &mut &[u8], state_count: usize) -> io::Result<(Vec<u8>, usize)> {
    let (context_size, is_compressed) = read_header(src)?;

    let len = read_uint7_as(src)?;
    let context_src = read_src(src, state_count, context_size, is_compressed)?;

    Ok((context_src, len))
}

fn read_header(src: &mut &[u8]) -> io::Result<(usize, bool)> {
    let n = read_uint7(src)?;

    let context_size =
        usize::try_from(n >> 1).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let is_compressed = (n & 0x01) == 0;

    Ok((context_size, is_compressed))
}

fn read_src(
    src: &mut &[u8],
    state_count: usize,
    len: usize,
    is_compressed: bool,
) -> io::Result<Vec<u8>> {
    if is_compressed {
        let comp_meta_len = read_uint7_as(src)?;

        let mut buf = vec![0; comp_meta_len];
        src.read_exact(&mut buf)?;

        let mut buf_reader = &buf[..];
        let mut dst = vec![0; len];
        order_0::decode(&mut buf_reader, &mut dst, state_count)?;

        Ok(dst)
    } else {
        let mut buf = vec![0; len];
        src.read_exact(&mut buf)?;
        Ok(buf)
    }
}

pub fn decode(mut src: &[u8], ctx: &[u8], len: usize) -> io::Result<Vec<u8>> {
    let mut context_src = ctx;

    let l = read_rle_table(&mut context_src)?;

    let mut dst = vec![0; len];
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
    let mut m = read_u8(src).map(u16::from)?;

    if m == 0 {
        m = 256;
    }

    let mut l = [false; 256];

    for _ in 0..m {
        let s = read_u8(src)?;
        l[usize::from(s)] = true;
    }

    Ok(l)
}
