use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let chunk_count = read_chunk_count(reader)?;

    let compressed_sizes = read_compressed_sizes(reader, chunk_count)?;
    let uncompressed_sizes = build_uncompressed_sizes(len, chunk_count);

    let chunks: Vec<_> = compressed_sizes
        .into_iter()
        .zip(uncompressed_sizes)
        .map(|(compressed_size, uncompressed_size)| {
            let mut buf = vec![0; compressed_size];
            reader.read_exact(&mut buf)?;
            super::decode(&mut &buf[..], uncompressed_size)
        })
        .collect::<io::Result<_>>()?;

    Ok(transpose(&chunks, len))
}

fn read_chunk_count<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    read_u8(reader).map(usize::from)
}

fn read_compressed_sizes<R>(reader: &mut R, chunk_count: usize) -> io::Result<Vec<usize>>
where
    R: Read,
{
    (0..chunk_count).map(|_| read_uint7_as(reader)).collect()
}

// § 3.6 "Striped rANS Nx16" (2023-03-15): "The uncompressed data length may not necessary be an
// exact multiple of _N_, in which case the latter uncompressed sub-streams may be 1 byte shorter."
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
        // § 3.6 "Striped rANS Nx16" (2023-03-15)
        let chunks = [&b"abcde"[..], &b"ABCD"[..], &b"1234"[..]];
        let actual = transpose(&chunks, 13);
        let expected = b"aA1bB2cC3dD4e";
        assert_eq!(actual, expected);
    }
}
