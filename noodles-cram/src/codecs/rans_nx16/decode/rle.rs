use std::io::{self, Cursor, Read};

use super::order_0;
use crate::io::reader::num::{read_u8, read_uint7};

pub(super) fn decode_rle_meta<R>(
    reader: &mut R,
    state_count: usize,
) -> io::Result<([bool; 256], Cursor<Vec<u8>>, usize)>
where
    R: Read,
{
    let rle_meta_len = read_uint7(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let len = read_uint7(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let rle_meta = if rle_meta_len & 1 == 1 {
        let mut buf = vec![0; rle_meta_len / 2];
        reader.read_exact(&mut buf)?;
        buf
    } else {
        let comp_meta_len = read_uint7(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let mut buf = vec![0; comp_meta_len];
        reader.read_exact(&mut buf)?;

        let mut buf_reader = &buf[..];
        let mut dst = vec![0; rle_meta_len / 2];
        order_0::decode(&mut buf_reader, &mut dst, state_count)?;

        dst
    };

    let mut rle_meta_reader = Cursor::new(rle_meta);

    let mut m = read_u8(&mut rle_meta_reader).map(u16::from)?;

    if m == 0 {
        m = 256;
    }

    let mut l = [false; 256];

    for _ in 0..m {
        let s = read_u8(&mut rle_meta_reader)?;
        l[usize::from(s)] = true;
    }

    Ok((l, rle_meta_reader, len))
}

pub fn decode<R>(
    mut src: &[u8],
    l: &[bool; 256],
    rle_meta: &mut R,
    len: usize,
) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut dst = vec![0; len];
    let mut j = 0;

    while j < dst.len() {
        let sym = read_u8(&mut src)?;

        if l[usize::from(sym)] {
            let run = read_uint7(rle_meta).and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

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
