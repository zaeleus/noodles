mod context;

use std::io;

use self::context::Context;
pub use self::context::read_context;

pub fn decode(src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    let mapping_table = ctx.mapping_table;

    let mut dst = vec![0; ctx.uncompressed_size];

    match ctx.symbol_count.get() {
        1 => dst.fill(mapping_table[0]),
        2 => unpack(src, mapping_table, 8, &mut dst),
        3..=4 => unpack(src, mapping_table, 4, &mut dst),
        5..=16 => unpack(src, mapping_table, 2, &mut dst),
        n => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected bit pack symbol count to be <= 16, got {n}"),
            ));
        }
    }

    Ok(dst)
}

fn unpack(src: &[u8], mapping_table: &[u8], chunk_size: usize, dst: &mut [u8]) {
    const BITS: usize = u8::BITS as usize;

    let shift = BITS / chunk_size;
    let mask = (1 << shift) - 1;

    for (mut s, chunk) in src.iter().copied().zip(dst.chunks_mut(chunk_size)) {
        for d in chunk {
            *d = mapping_table[usize::from(s & mask)];
            s >>= shift;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unpack() {
        let mut dst = [0; 4];
        unpack(&[0b00001100], b"nd", 8, &mut dst);
        assert_eq!(&dst, b"nndd");

        let mut dst = [0; 4];
        unpack(&[0b11100100], b"ndls", 4, &mut dst);
        assert_eq!(&dst, b"ndls");

        let mut dst = [0; 2];
        unpack(&[0b00100000], b"nodles", 2, &mut dst);
        assert_eq!(&dst, b"nd");
    }
}
