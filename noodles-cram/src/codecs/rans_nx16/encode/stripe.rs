use std::io::{self, Write};

use crate::{
    codecs::rans_nx16::Flags,
    io::writer::num::{write_u8, write_uint7},
};

const CHUNK_COUNT: usize = 4;

pub(super) fn rans_encode_stripe(src: &[u8]) -> io::Result<Vec<u8>> {
    let mut ulens = Vec::with_capacity(CHUNK_COUNT);

    for j in 0..CHUNK_COUNT {
        let mut ulen = src.len() / CHUNK_COUNT;

        if src.len() % CHUNK_COUNT > j {
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
    write_chunk_count(&mut dst, compressed_chunks.len())?;
    write_compressed_sizes(&mut dst, &compressed_chunks)?;

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

fn write_chunk_count(dst: &mut Vec<u8>, chunk_count: usize) -> io::Result<()> {
    let n =
        u8::try_from(chunk_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

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
        write_chunk_count(&mut dst, 4)?;
        assert_eq!(dst, [0x04]);

        dst.clear();
        assert!(matches!(
            write_chunk_count(&mut dst, 256),
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
