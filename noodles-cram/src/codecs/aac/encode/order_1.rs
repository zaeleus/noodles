use std::{io, num::NonZero};

use super::write_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    let symbol_count = NonZero::new(usize::from(max_sym) + 1).unwrap();
    write_symbol_count(dst, symbol_count)?;

    let mut models = vec![Model::new(symbol_count); symbol_count.get()];

    let mut range_coder = RangeCoder::default();

    models[0].encode(dst, &mut range_coder, src[0])?;

    for window in src.windows(2) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = window[1];
        models[sym_0].encode(dst, &mut range_coder, sym_1)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
