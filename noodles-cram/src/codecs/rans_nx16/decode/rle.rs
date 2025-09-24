use std::io::{self, Read};

use super::order_0;
use crate::io::reader::num::{read_u8, read_uint7, read_uint7_as};

pub(super) fn decode_rle_meta<R>(reader: &mut R, state_count: usize) -> io::Result<(Vec<u8>, usize)>
where
    R: Read,
{
    let (context_size, is_compressed) = read_header(reader)?;

    let len = read_uint7_as(reader)?;
    let context_src = read_src(reader, state_count, context_size, is_compressed)?;

    Ok((context_src, len))
}

fn read_header<R>(reader: &mut R) -> io::Result<(usize, bool)>
where
    R: Read,
{
    let n = read_uint7(reader)?;

    let context_size =
        usize::try_from(n >> 1).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let is_compressed = (n & 0x01) == 0;

    Ok((context_size, is_compressed))
}

fn read_src<R>(
    reader: &mut R,
    state_count: usize,
    len: usize,
    is_compressed: bool,
) -> io::Result<Vec<u8>>
where
    R: Read,
{
    if is_compressed {
        let comp_meta_len = read_uint7_as(reader)?;

        let mut buf = vec![0; comp_meta_len];
        reader.read_exact(&mut buf)?;

        let mut buf_reader = &buf[..];
        let mut dst = vec![0; len];
        order_0::decode(&mut buf_reader, &mut dst, state_count)?;

        Ok(dst)
    } else {
        let mut buf = vec![0; len];
        reader.read_exact(&mut buf)?;
        Ok(buf)
    }
}

pub fn decode<R>(mut src: &[u8], rle_meta: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let l = read_rle_table(rle_meta)?;

    let mut dst = vec![0; len];
    let mut j = 0;

    while j < dst.len() {
        let sym = read_u8(&mut src)?;

        if l[usize::from(sym)] {
            let run = read_uint7_as(rle_meta)?;

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

fn read_rle_table<R>(reader: &mut R) -> io::Result<[bool; 256]>
where
    R: Read,
{
    let mut m = read_u8(reader).map(u16::from)?;

    if m == 0 {
        m = 256;
    }

    let mut l = [false; 256];

    for _ in 0..m {
        let s = read_u8(reader)?;
        l[usize::from(s)] = true;
    }

    Ok(l)
}
