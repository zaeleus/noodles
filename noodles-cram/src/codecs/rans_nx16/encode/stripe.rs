use std::io::{self, Write};

use crate::{
    codecs::rans_nx16::Flags,
    io::writer::num::{write_u8, write_uint7},
};

pub(super) fn rans_encode_stripe(src: &[u8], n: usize) -> io::Result<Vec<u8>> {
    let mut ulens = Vec::with_capacity(n);

    for j in 0..n {
        let mut ulen = src.len() / n;

        if src.len() % n > j {
            ulen += 1;
        }

        ulens.push(ulen);
    }

    let uncompressed_chunks = transpose(src, &ulens);

    let compressed_chunks: Vec<_> = uncompressed_chunks
        .into_iter()
        .map(|chunk| super::encode(Flags::empty(), &chunk))
        .collect::<io::Result<_>>()?;

    let mut dst = Vec::new();

    write_u8(&mut dst, n as u8)?;

    for chunk in &compressed_chunks {
        let clen = chunk.len() as u32;
        write_uint7(&mut dst, clen)?;
    }

    for chunk in &compressed_chunks {
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
