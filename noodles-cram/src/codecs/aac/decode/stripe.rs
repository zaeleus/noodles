use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode(src: &mut &[u8], len: usize) -> io::Result<Vec<u8>> {
    let chunk_count = read_chunk_count(src)?;

    let compressed_sizes = read_compressed_sizes(src, chunk_count)?;
    let uncompressed_sizes = build_uncompressed_sizes(len, chunk_count);

    let mut t = Vec::with_capacity(chunk_count);

    for (compressed_size, uncompressed_size) in compressed_sizes.into_iter().zip(uncompressed_sizes)
    {
        let mut buf = vec![0; compressed_size];
        src.read_exact(&mut buf)?;
        let chunk = super::decode(&buf, uncompressed_size)?;
        t.push(chunk);
    }

    Ok(transpose(&t, len))
}

fn read_chunk_count(src: &mut &[u8]) -> io::Result<usize> {
    read_u8(src).map(usize::from)
}

fn read_compressed_sizes(src: &mut &[u8], chunk_count: usize) -> io::Result<Vec<usize>> {
    (0..chunk_count).map(|_| read_uint7_as(src)).collect()
}

fn build_uncompressed_sizes(len: usize, chunk_count: usize) -> Vec<usize> {
    (0..chunk_count)
        .map(|i| {
            let (q, r) = (len / chunk_count, len % chunk_count);
            if r > i { q + 1 } else { q }
        })
        .collect()
}

fn transpose<T>(chunks: &[T], uncompressed_size: usize) -> Vec<u8>
where
    T: AsRef<[u8]>,
{
    let mut dst = vec![0; uncompressed_size];

    for (i, chunk) in chunks.iter().enumerate() {
        for (j, s) in chunk.as_ref().iter().enumerate() {
            dst[j * chunks.len() + i] = *s;
        }
    }

    dst
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_uncompressed_sizes() {
        assert_eq!(build_uncompressed_sizes(13, 3), [5, 4, 4]);
    }

    #[test]
    fn test_transpose() {
        let chunks = [&b"abcde"[..], &b"ABCD"[..], &b"1234"[..]];
        let actual = transpose(&chunks, 13);
        let expected = b"aA1bB2cC3dD4e";
        assert_eq!(actual, expected);
    }
}
