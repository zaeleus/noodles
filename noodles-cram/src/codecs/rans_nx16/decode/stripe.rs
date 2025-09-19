use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let x = read_u8(reader).map(usize::from)?;
    let mut clens: Vec<usize> = Vec::with_capacity(x);

    for _ in 0..x {
        let clen = read_uint7_as(reader)?;
        clens.push(clen);
    }

    let mut ulens = Vec::with_capacity(x);
    let mut t = Vec::with_capacity(x);

    for j in 0..x {
        let mut ulen = len / x;

        if len % x > j {
            ulen += 1;
        }

        let chunk = super::decode(reader, ulen)?;

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut dst = vec![0; len];

    for j in 0..x {
        for i in 0..ulens[j] {
            dst[i * x + j] = t[j][i];
        }
    }

    Ok(dst)
}
