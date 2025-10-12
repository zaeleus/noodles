use std::io::{self, Write};

use crate::{
    codecs::rans_nx16::Flags,
    io::writer::num::{write_u8, write_uint7},
};

pub(super) fn rans_encode_stripe(src: &[u8], n: usize) -> io::Result<Vec<u8>> {
    let mut ulens = Vec::with_capacity(n);
    let mut t = Vec::with_capacity(n);

    for j in 0..n {
        let mut ulen = src.len() / n;

        if src.len() % n > j {
            ulen += 1;
        }

        let chunk = vec![0; ulen];

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut x = 0;
    let mut i = 0;

    while i < src.len() {
        for j in 0..n {
            if x < ulens[j] {
                t[j][x] = src[i + j];
            }
        }

        x += 1;
        i += n;
    }

    let mut chunks = vec![Vec::new(); n];

    for (chunk, s) in chunks.iter_mut().zip(t.iter()) {
        *chunk = super::encode(Flags::empty(), s)?;
    }

    let mut dst = Vec::new();

    write_u8(&mut dst, n as u8)?;

    for chunk in &chunks {
        let clen = chunk.len() as u32;
        write_uint7(&mut dst, clen)?;
    }

    for chunk in &chunks {
        dst.write_all(chunk)?;
    }

    Ok(dst)
}
