use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode(src: &mut &[u8], len: usize) -> io::Result<Vec<u8>> {
    let n = read_u8(src).map(usize::from)?;
    let mut clens: Vec<usize> = Vec::with_capacity(n);

    for _ in 0..n {
        let clen = read_uint7_as(src)?;
        clens.push(clen);
    }

    let mut ulens = Vec::with_capacity(n);
    let mut t = Vec::with_capacity(n);

    for (j, clen) in clens.iter().enumerate() {
        let mut ulen = len / n;

        if len % n > j {
            ulen += 1;
        }

        let mut buf = vec![0; *clen];
        src.read_exact(&mut buf)?;
        let chunk = super::decode(&buf, ulen)?;

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut dst = vec![0; len];

    for j in 0..n {
        for i in 0..ulens[j] {
            dst[i * n + j] = t[j][i];
        }
    }

    Ok(dst)
}
