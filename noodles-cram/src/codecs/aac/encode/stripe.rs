use std::io::{self, Write};

use crate::{
    codecs::aac::Flags,
    io::writer::num::{write_u8, write_uint7},
};

pub(super) fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    const N: usize = 4;

    let mut ulens = Vec::with_capacity(N);

    for j in 0..N {
        let mut ulen = src.len() / N;

        if src.len() % N > j {
            ulen += 1;
        }

        ulens.push(ulen);
    }

    let mut chunks = vec![Vec::new(); N];
    let t = transpose(src, &ulens);

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

fn transpose(src: &[u8], chunk_sizes: &[usize]) -> Vec<Vec<u8>> {
    let mut chunks: Vec<_> = chunk_sizes.iter().map(|&len| vec![0; len]).collect();

    for (i, chunk) in chunks.iter_mut().enumerate() {
        for (j, d) in chunk.iter_mut().enumerate() {
            *d = src[j * chunk_sizes.len() + i];
        }
    }

    chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transpose() {
        let chunk_sizes = [5, 4, 4];
        let src = b"aA1bB2cC3dD4e";

        assert_eq!(
            transpose(src, &chunk_sizes),
            [Vec::from(b"abcde"), Vec::from(b"ABCD"), Vec::from(b"1234")]
        )
    }
}
