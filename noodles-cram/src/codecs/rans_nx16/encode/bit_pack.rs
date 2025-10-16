pub mod context;

pub use self::context::{Context, build_context, write_context};
use crate::codecs::rans_nx16::ALPHABET_SIZE;

pub fn encode(src: &[u8], ctx: &Context) -> Vec<u8> {
    let mapping_table = &ctx.mapping_table;

    match ctx.symbol_count.get() {
        1 => Vec::new(),
        2 => pack(src, mapping_table, 8),
        3..=4 => pack(src, mapping_table, 4),
        5..=16 => pack(src, mapping_table, 2),
        _ => unreachable!(),
    }
}

fn pack(src: &[u8], mapping_table: &[u8; ALPHABET_SIZE], chunk_size: usize) -> Vec<u8> {
    const BITS: usize = u8::BITS as usize;

    let shift = BITS / chunk_size;

    let compressed_size = src.len().div_ceil(chunk_size);
    let mut dst = vec![0; compressed_size];

    for (d, chunk) in dst.iter_mut().zip(src.chunks(chunk_size)) {
        for (i, &sym) in chunk.iter().enumerate() {
            let n = mapping_table[usize::from(sym)];
            *d |= n << (shift * i);
        }
    }

    dst
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_8() {
        let src = [b'n', b'n', b'd', b'd'];

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'd')] = 1;

        assert_eq!(pack(&src, &mapping_table, 8), [0b00001100]);
    }

    #[test]
    fn test_pack_4() {
        let src = [b'n', b'd', b'l', b's'];

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'd')] = 1;
        mapping_table[usize::from(b'l')] = 2;
        mapping_table[usize::from(b's')] = 3;

        assert_eq!(pack(&src, &mapping_table, 4), [0b11100100]);
    }

    #[test]
    fn test_pack_2() {
        let src = [b'n', b'd'];

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'o')] = 1;
        mapping_table[usize::from(b'd')] = 2;
        mapping_table[usize::from(b'l')] = 3;
        mapping_table[usize::from(b's')] = 4;

        assert_eq!(pack(&src, &mapping_table, 2), [0b00100000]);
    }
}
