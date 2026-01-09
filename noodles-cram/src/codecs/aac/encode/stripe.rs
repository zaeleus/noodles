use std::io::{self, Write};

use crate::{
    codecs::aac::Flags,
    io::writer::num::{write_u8, write_uint7},
};

pub(super) fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    const N: usize = 4;

    let mut ulens = Vec::with_capacity(N);
    let mut t = Vec::with_capacity(N);

    for j in 0..N {
        let mut ulen = src.len() / N;

        if src.len() % N > j {
            ulen += 1;
        }

        let chunk = vec![0; ulen];

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut x = 0;
    let mut i = 0;

    while i < src.len() {
        for j in 0..N {
            if x < ulens[j] {
                t[j][x] = src[i + j];
            }
        }

        x += 1;
        i += N;
    }

    let mut chunks = vec![Vec::new(); N];

    for (chunk, s) in chunks.iter_mut().zip(t.iter()) {
        *chunk = super::encode(Flags::empty(), s)?;
    }

    let mut dst = Vec::new();

    write_u8(&mut dst, N as u8)?;

    for chunk in &chunks {
        let clen = chunk.len() as u32;
        write_uint7(&mut dst, clen)?;
    }

    for chunk in &chunks {
        dst.write_all(chunk)?;
    }

    Ok(dst)
}
