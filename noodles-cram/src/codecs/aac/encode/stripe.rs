use std::{io, num::NonZero};

use crate::{
    codecs::aac::Flags,
    io::writer::num::{write_u8, write_uint7},
};

const CHUNK_COUNT: NonZero<usize> = NonZero::new(4).unwrap();

pub(super) fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    let mut ulens = Vec::with_capacity(CHUNK_COUNT.get());

    for j in 0..CHUNK_COUNT.get() {
        // SAFETY: `CHUNK_COUNT` > 0.
        let mut ulen = src.len() / CHUNK_COUNT.get();

        // SAFETY: `CHUNK_COUNT` > 0.
        if src.len() % CHUNK_COUNT.get() > j {
            ulen += 1;
        }

        ulens.push(ulen);
    }

    let mut chunks = vec![Vec::new(); CHUNK_COUNT.get()];
    let t = transpose(src, &ulens);

    for (chunk, s) in chunks.iter_mut().zip(t.iter()) {
        *chunk = super::encode(Flags::NO_SIZE, s)?;
    }

    let mut dst = Vec::new();
    write_chunk_count(&mut dst, CHUNK_COUNT)?;
    write_compressed_sizes(&mut dst, &chunks)?;
    dst.extend(chunks.into_iter().flatten());

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
