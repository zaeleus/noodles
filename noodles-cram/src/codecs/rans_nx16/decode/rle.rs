use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_uint7};

pub fn decode<R>(
    mut src: &[u8],
    l: &[bool; 256],
    rle_meta: &mut R,
    len: usize,
) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut dst = vec![0; len];
    let mut j = 0;

    while j < dst.len() {
        let sym = read_u8(&mut src)?;

        if l[usize::from(sym)] {
            let run = read_uint7(rle_meta).and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

            for k in 0..=run {
                dst[j + k] = sym;
            }

            j += run + 1;
        } else {
            dst[j] = sym;
            j += 1;
        }
    }

    Ok(dst)
}
