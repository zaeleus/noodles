use std::io;

use super::{count_symbols, write_symbol_count};
use crate::codecs::aac::{Model, RangeCoder};

const CONTEXT_SIZE: usize = 2;
const NUL: u8 = 0x00;

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = count_symbols(src);
    write_symbol_count(dst, symbol_count)?;

    let mut models = vec![Model::new(symbol_count); symbol_count.get()];
    let mut coder = RangeCoder::default();

    models[usize::from(NUL)].encode(dst, &mut coder, src[0])?;

    for syms in src.windows(CONTEXT_SIZE) {
        let (prev_sym, sym) = (usize::from(syms[0]), syms[1]);
        models[prev_sym].encode(dst, &mut coder, sym)?;
    }

    coder.range_encode_end(dst)?;

    Ok(())
}
