use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let chunk_count = read_chunk_count(reader)?;

    let _compressed_sizes = read_compressed_sizes(reader, chunk_count)?;

    let mut t = Vec::with_capacity(chunk_count);

    for j in 0..chunk_count {
        let mut ulen = len / chunk_count;

        if len % chunk_count > j {
            ulen += 1;
        }

        let chunk = super::decode(reader, ulen)?;

        t.push(chunk);
    }

    Ok(transpose(&t, len))
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
    fn test_transpose() {
        // § 3.6 "Striped rANS Nx16" (2023-03-15)
        let chunks = [&b"abcde"[..], &b"ABCD"[..], &b"1234"[..]];
        let actual = transpose(&chunks, 13);
        let expected = b"aA1bB2cC3dD4e";
        assert_eq!(actual, expected);
    }
}
