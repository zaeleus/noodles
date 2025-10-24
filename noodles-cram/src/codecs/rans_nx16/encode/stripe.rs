use std::{io, num::NonZero};

use crate::{
    codecs::rans_nx16::Flags,
    io::writer::num::{write_u8, write_uint7},
};

const CHUNK_COUNT: NonZero<usize> = NonZero::new(4).unwrap();

pub(super) fn rans_encode_stripe(src: &[u8]) -> io::Result<Vec<u8>> {
    let uncompressed_sizes = build_uncompressed_sizes(src.len(), CHUNK_COUNT);
    let uncompressed_chunks = transpose(src, &uncompressed_sizes);

    let compressed_chunks: Vec<_> = uncompressed_chunks
        .into_iter()
        .map(|chunk| super::encode(Flags::empty(), &chunk))
        .collect::<io::Result<_>>()?;

    let mut dst = Vec::new();
    write_chunk_count(&mut dst, CHUNK_COUNT)?;
    write_compressed_sizes(&mut dst, &compressed_chunks)?;
    dst.extend(compressed_chunks.into_iter().flatten());

    Ok(dst)
}

// § 3.6 "Striped rANS Nx16" (2023-03-15): "The uncompressed data length may not necessary be an
// exact multiple of _N_, in which case the latter uncompressed sub-streams may be 1 byte shorter."
fn build_uncompressed_sizes(len: usize, chunk_count: NonZero<usize>) -> Vec<usize> {
    let n = chunk_count.get();

    (0..n)
        .map(|i| {
            // SAFETY: `n` > 0.
            let (q, r) = (len / n, len % n);
            if r > i { q + 1 } else { q }
        })
        .collect()
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

fn write_chunk_count(dst: &mut Vec<u8>, chunk_count: NonZero<usize>) -> io::Result<()> {
    let n = u8::try_from(chunk_count.get())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_u8(dst, n)?;

    Ok(())
}

fn write_compressed_sizes(dst: &mut Vec<u8>, chunks: &[Vec<u8>]) -> io::Result<()> {
    for chunk in chunks {
        let n = u32::try_from(chunk.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        write_uint7(dst, n)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_uncompressed_sizes() {
        assert_eq!(
            build_uncompressed_sizes(13, const { NonZero::new(3).unwrap() }),
            [5, 4, 4]
        );
    }

    #[test]
    fn test_transpose() {
        let chunk_sizes = [5, 4, 4];
        let src = b"aA1bB2cC3dD4e";

        assert_eq!(
            transpose(src, &chunk_sizes),
            [Vec::from(b"abcde"), Vec::from(b"ABCD"), Vec::from(b"1234")]
        )
    }

    #[test]
    fn test_write_chunk_count() -> io::Result<()> {
        let mut dst = Vec::new();

        dst.clear();
        write_chunk_count(&mut dst, const { NonZero::new(4).unwrap() })?;
        assert_eq!(dst, [0x04]);

        dst.clear();
        assert!(matches!(
            write_chunk_count(&mut dst, const { NonZero::new(256).unwrap() }),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_compressed_sizes() -> io::Result<()> {
        let mut dst = Vec::new();
        let compressed_chunks = [Vec::from(b"abcde"), Vec::from(b"ABCD"), Vec::from(b"1234")];
        write_compressed_sizes(&mut dst, &compressed_chunks)?;
        assert_eq!(dst, [0x05, 0x04, 0x04]);
        Ok(())
    }
}
