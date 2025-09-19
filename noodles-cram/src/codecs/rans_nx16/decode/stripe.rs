use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub(super) fn decode<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let x = read_u8(reader).map(usize::from)?;
    let mut clens: Vec<usize> = Vec::with_capacity(x);

    for _ in 0..x {
        let clen = read_uint7_as(reader)?;
        clens.push(clen);
    }

    let mut t = Vec::with_capacity(x);

    for j in 0..x {
        let mut ulen = len / x;

        if len % x > j {
            ulen += 1;
        }

        let chunk = super::decode(reader, ulen)?;

        t.push(chunk);
    }

    Ok(transpose(&t, len))
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
