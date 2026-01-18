use std::io;

use super::{count_symbols, write_symbol_count};
use crate::codecs::aac::{Model, RangeCoder};

const CONTEXT_SIZE: usize = 2;

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = count_symbols(src);
    write_symbol_count(dst, symbol_count)?;

    let mut models = vec![Model::new(symbol_count); symbol_count.get()];

    let mut range_coder = RangeCoder::default();

    models[0].encode(dst, &mut range_coder, src[0])?;

    for window in src.windows(CONTEXT_SIZE) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = window[1];
        models[sym_0].encode(dst, &mut range_coder, sym_1)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
