mod context;

use std::io;

pub use self::context::{Context, read_context};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::reader::num::{read_u8, read_uint7_as},
};

pub fn decode(mut src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    let mut context_src = &ctx.src[..];

    let rle_alphabet = read_rle_alphabet(&mut context_src)?;

    let mut dst = vec![0; ctx.len];
    let mut iter = dst.iter_mut();

    while let Some(d) = iter.next() {
        let sym = read_u8(&mut src)?;

        *d = sym;

        if rle_alphabet[usize::from(sym)] {
            let len = read_uint7_as(&mut context_src)?;

            for e in iter.by_ref().take(len) {
                *e = sym;
            }
        }
    }

    Ok(dst)
}

fn read_rle_alphabet(src: &mut &[u8]) -> io::Result<[bool; ALPHABET_SIZE]> {
    let mut alphabet = [false; ALPHABET_SIZE];

    let symbol_count = {
        let n = read_u8(src).map(usize::from)?;
        if n == 0 { ALPHABET_SIZE } else { n }
    };

    for _ in 0..symbol_count {
        let sym = read_u8(src)?;
        alphabet[usize::from(sym)] = true;
    }

    Ok(alphabet)
}
