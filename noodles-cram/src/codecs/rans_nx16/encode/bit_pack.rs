mod context;

use std::io;

pub use self::context::{Context, build_context};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::writer::num::{write_u8, write_uint7},
};

pub fn encode(src: &[u8], ctx: &Context) -> io::Result<(Vec<u8>, Vec<u8>)> {
    let mapping_table = &ctx.mapping_table;

    let buf = match ctx.symbol_count.get() {
        1 => Vec::new(),
        2 => pack(src, mapping_table, 8),
        3..=4 => pack(src, mapping_table, 4),
        5..=16 => pack(src, mapping_table, 2),
        _ => unreachable!(),
    };

    let mut header = Vec::new();

    let n = u8::try_from(ctx.symbol_count.get())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u8(&mut header, n)?;

    for (sym, _) in ctx.alphabet.iter().enumerate().filter(|(_, a)| **a) {
        write_u8(&mut header, sym as u8)?;
    }

    let len = buf.len() as u32;
    write_uint7(&mut header, len)?;

    Ok((header, buf))
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
