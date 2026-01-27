use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode(src: &mut &[u8], len: usize) -> io::Result<Vec<u8>> {
    let chunk_count = read_chunk_count(src)?;
    let mut clens: Vec<usize> = Vec::with_capacity(chunk_count);

    for _ in 0..chunk_count {
        let clen = read_uint7_as(src)?;
        clens.push(clen);
    }

    let mut t = Vec::with_capacity(chunk_count);

    for (j, clen) in clens.iter().enumerate() {
        let mut ulen = len / chunk_count;

        if len % chunk_count > j {
            ulen += 1;
        }

        let mut buf = vec![0; *clen];
        src.read_exact(&mut buf)?;
        let chunk = super::decode(&buf, ulen)?;

        t.push(chunk);
    }

    Ok(transpose(&t, len))
}

fn read_chunk_count(src: &mut &[u8]) -> io::Result<usize> {
    read_u8(src).map(usize::from)
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
        let chunks = [&b"abcde"[..], &b"ABCD"[..], &b"1234"[..]];
        let actual = transpose(&chunks, 13);
        let expected = b"aA1bB2cC3dD4e";
        assert_eq!(actual, expected);
    }
}
